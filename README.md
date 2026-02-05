# CardioTwin (Combined Final App)

This workspace now combines both hackathon worktrees into a single runnable app:

- `backend/`: FastAPI API backed by converted MATLAB cardiac simulation kernels.
- `frontend/`: React + Vite interactive simulation UI.
- `models/`: source MATLAB assets and calibration files consumed by the backend.

## Run Backend

```bash
cd /Users/dev/Documents/hackathons/cerebral_valley/codex/backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

Backend endpoints:

- `GET /health`
- `GET /state`
- `POST /simulate/step`
- `POST /interventions`
- `POST /reset`

## Run Frontend

```bash
cd /Users/dev/Documents/hackathons/cerebral_valley/codex/frontend
npm install
npm run dev
```

Optional API base override:

```bash
VITE_API_BASE_URL=http://localhost:8000 npm run dev
```

## Notes

- The frontend contract remains stable while the backend now runs converted model dynamics.
- Intervention controls in the UI map to converted model parameters through an adapter layer in `backend/app/simulation/engine.py`.
