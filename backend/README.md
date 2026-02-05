# CardioTwin Backend (Converted Model)

FastAPI service exposing the MVP contract from `hackthon_approach.md`, now backed by converted MATLAB model kernels.

## Endpoints
- `GET /state`
- `POST /simulate/step`
- `POST /interventions`
- `POST /reset`
- `GET /health`

## Run
```bash
cd backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

The backend simulation engine is an adapter over the converted MATLAB ODE system with the frontend schema preserved.
