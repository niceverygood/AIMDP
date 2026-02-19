"""
Training script for property prediction models.
Trains a neural network head on top of Uni-Mol embeddings.

Usage:
    python -m ml.train --dataset qm9 --epochs 50 --batch-size 64
"""

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset, random_split

logger = logging.getLogger(__name__)

# Properties to predict
PROPERTY_NAMES = [
    "thermal_stability",
    "dielectric_constant",
    "bandgap",
    "solubility",
    "density",
]


class PropertyPredictor(nn.Module):
    """
    Neural network head for property prediction.
    Input: Uni-Mol embedding (768-dim)
    Output: 5 property values
    """

    def __init__(self, input_dim: int = 768, num_properties: int = 5):
        super().__init__()
        self.head = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 128),
            nn.LayerNorm(128),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, num_properties),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(x)


class ConditionalMolVAE(nn.Module):
    """
    Conditional VAE for inverse molecular design.
    Given target properties, generate SMILES of candidate molecules.
    """

    def __init__(
        self,
        vocab_size: int = 128,
        latent_dim: int = 128,
        condition_dim: int = 5,
        hidden_dim: int = 256,
    ):
        super().__init__()
        self.latent_dim = latent_dim

        # Encoder: SMILES tokens → latent space
        self.encoder = nn.GRU(
            vocab_size, hidden_dim, batch_first=True, bidirectional=True
        )
        self.fc_mu = nn.Linear(hidden_dim * 2, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim * 2, latent_dim)

        # Decoder: latent + condition → SMILES tokens
        self.decoder = nn.GRU(
            latent_dim + condition_dim, hidden_dim, batch_first=True
        )
        self.output = nn.Linear(hidden_dim, vocab_size)

    def encode(self, x: torch.Tensor):
        _, h = self.encoder(x)
        h = torch.cat([h[0], h[1]], dim=-1)  # Combine bidirectional
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(
        self, z: torch.Tensor, condition: torch.Tensor, seq_len: int = 100
    ) -> torch.Tensor:
        combined = torch.cat([z, condition], dim=-1).unsqueeze(1).repeat(1, seq_len, 1)
        output, _ = self.decoder(combined)
        return self.output(output)

    def forward(self, x: torch.Tensor, condition: torch.Tensor):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z, condition, x.size(1))
        return recon, mu, logvar

    def generate(self, target_properties: list[float], num_samples: int = 10):
        """Generate molecules conditioned on target properties."""
        self.eval()
        with torch.no_grad():
            z = torch.randn(num_samples, self.latent_dim)
            condition = torch.tensor(target_properties).unsqueeze(0).repeat(num_samples, 1)
            output = self.decode(z, condition)
            tokens = output.argmax(dim=-1)
        return tokens


def train_property_predictor(
    embeddings: np.ndarray,
    properties: np.ndarray,
    epochs: int = 50,
    batch_size: int = 64,
    learning_rate: float = 1e-3,
    save_path: str = "./ml/property_models/predictor.pt",
):
    """
    Train the property prediction model.

    Args:
        embeddings: (N, 768) array of molecular embeddings
        properties: (N, 5) array of property values
        epochs: Number of training epochs
        batch_size: Batch size
        learning_rate: Learning rate
        save_path: Path to save the trained model
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Training on {device}")

    # Prepare data
    X = torch.tensor(embeddings, dtype=torch.float32)
    y = torch.tensor(properties, dtype=torch.float32)

    dataset = TensorDataset(X, y)
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = random_split(dataset, [train_size, val_size])

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)

    # Model
    model = PropertyPredictor(input_dim=embeddings.shape[1]).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    best_val_loss = float("inf")

    for epoch in range(epochs):
        # Train
        model.train()
        train_loss = 0
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            optimizer.zero_grad()
            output = model(batch_X)
            loss = criterion(output, batch_y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        # Validate
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch_X, batch_y in val_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                output = model(batch_X)
                loss = criterion(output, batch_y)
                val_loss += loss.item()

        train_loss /= len(train_loader)
        val_loss /= len(val_loader)
        scheduler.step()

        if (epoch + 1) % 5 == 0:
            logger.info(
                f"Epoch {epoch+1}/{epochs} — "
                f"Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}"
            )

        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            torch.save(model, save_path)
            logger.info(f"Saved best model (val_loss={val_loss:.4f})")

    logger.info(f"Training complete. Best val loss: {best_val_loss:.4f}")
    return model


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train property prediction model")
    parser.add_argument("--dataset", default="qm9", help="Dataset to use")
    parser.add_argument("--epochs", type=int, default=50, help="Number of epochs")
    parser.add_argument("--batch-size", type=int, default=64, help="Batch size")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--save-path", default="./ml/property_models/predictor.pt")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # For demo: generate synthetic data
    logger.info("Generating synthetic training data for demo...")
    N = 1000
    embeddings = np.random.randn(N, 768).astype(np.float32)
    properties = np.random.randn(N, 5).astype(np.float32)

    train_property_predictor(
        embeddings=embeddings,
        properties=properties,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        save_path=args.save_path,
    )
