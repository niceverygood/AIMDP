# Materials AI Platform

AI 기반 신소재 탐색 플랫폼 — AI-powered Materials Discovery Platform

연구원이 원하는 물성을 입력하면 AI가 후보 소재를 추천하고, 데이터 파이프라인을 통해 원시 데이터를 학습 가능한 데이터셋으로 변환합니다.

## 핵심 기능

### 1. AI 소재 탐색 엔진
- 물성 조건 입력 (열안정성, 유전율, 밴드갭, 용해도 등)
- 벡터 유사도 검색 (pgvector + Uni-Mol 768차원 임베딩)
- 다목적 최적화 랭킹으로 최적 후보 추천
- 3D 분자 시각화 + 물성 레이더 차트

### 2. 데이터 파이프라인
- 수집 → 정규화 → 정제 → 특성 추출 → 임베딩 → 데이터셋
- WebSocket 기반 실시간 진행률 모니터링
- CSV, SDF, MOL, JSON 다양한 입력 형식 지원

### 3. 실험 피드백 루프 (Active Learning)
- 연구원 실험 결과 → 모델 재학습
- 예측값 vs 실측값 오차 분석
- 자동 재학습 트리거

## 기술 스택

| 영역 | 기술 |
|------|------|
| **Frontend** | Next.js 14 (App Router) + TypeScript + Tailwind CSS + shadcn/ui |
| **Backend** | FastAPI (Python 3.11+) |
| **AI/ML** | PyTorch + Uni-Mol + RDKit |
| **Database** | PostgreSQL 16 + pgvector |
| **Cache** | Redis 7 |
| **Infra** | Docker Compose |

## 빠른 시작

### 사전 요구사항
- Docker & Docker Compose
- Node.js 20+ (프론트엔드 개발시)
- Python 3.11+ (백엔드 개발시)

### Docker Compose로 실행

```bash
# 1. 환경변수 설정
cp .env.example .env

# 2. 전체 서비스 실행
docker compose up -d

# 3. 접속
# Frontend: http://localhost:3000
# Backend API: http://localhost:8000
# API Docs: http://localhost:8000/docs
```

### 로컬 개발 환경

```bash
# Backend
cd backend
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000

# Frontend (별도 터미널)
cd frontend
npm install
npm run dev
```

## 프로젝트 구조

```
materials-ai/
├── frontend/                     # Next.js 14 앱
│   ├── app/
│   │   ├── page.tsx              # 대시보드
│   │   ├── search/page.tsx       # AI 소재 탐색
│   │   └── pipeline/page.tsx     # 데이터 파이프라인
│   ├── components/
│   │   ├── ui/                   # shadcn 컴포넌트
│   │   ├── MoleculeViewer.tsx    # 3D 분자 뷰어
│   │   ├── PropertyRadar.tsx     # 물성 레이더 차트
│   │   ├── SearchFilters.tsx     # 탐색 조건 패널
│   │   ├── ResultCard.tsx        # 결과 카드
│   │   └── PipelineStep.tsx      # 파이프라인 단계
│   └── lib/
│       ├── api.ts                # API 클라이언트
│       └── types.ts              # TypeScript 타입
├── backend/
│   ├── app/
│   │   ├── main.py               # FastAPI 진입점
│   │   ├── routers/              # API 라우터
│   │   ├── services/             # 비즈니스 로직
│   │   ├── models/               # DB 모델 + 스키마
│   │   └── core/                 # 설정 + DB 연결
│   ├── ml/
│   │   ├── unimol/               # Uni-Mol 래퍼
│   │   ├── train.py              # 학습 스크립트
│   │   └── inference.py          # 추론 파이프라인
│   └── data/
│       ├── processors/           # 데이터 전처리
│       └── loaders/              # 파일 로더
├── docker-compose.yml
└── .env.example
```

## API 엔드포인트

| Method | Endpoint | 설명 |
|--------|----------|------|
| `POST` | `/api/search` | AI 소재 탐색 |
| `POST` | `/api/search/predict` | 물성 예측 |
| `GET` | `/api/materials` | 소재 목록 |
| `GET` | `/api/materials/{id}` | 소재 상세 |
| `POST` | `/api/materials` | 소재 등록 |
| `POST` | `/api/materials/feedback` | 실험 피드백 |
| `GET` | `/api/materials/stats` | DB 통계 |
| `POST` | `/api/pipeline/start` | 파이프라인 시작 |
| `GET` | `/api/pipeline/{id}` | 파이프라인 상태 |
| `WS` | `/api/pipeline/ws/{id}` | 실시간 진행률 |

## 디자인 시스템

- **테마**: 다크 네이비 (#0A1628, #0D1B2A)
- **액센트**: 시안 (#00B4D8), 그린 (#00E676), 퍼플 (#7C4DFF)
- **차트**: Recharts + 3Dmol.js
- **애니메이션**: Framer Motion
