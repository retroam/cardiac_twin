# CardioTwin

CardioTwin is an interactive cardiac simulation app built for the OpenAI Codex Hackathon.
It lets users inhibit signaling pathways and observe synchronized downstream effects across:

- molecular signaling activity
- electrical behavior (ECG-like trace)
- mechanical pumping behavior

The project combines two parallel build streams:

- converted MATLAB model kernels for backend dynamics
- interactive React simulation UI for intervention and visualization

## Project Goals

- Demonstrate causal pathway intervention, not static visualization.
- Show real-time propagation from signaling changes to electrical and mechanical outputs.
- Keep a stable frontend/backend API contract for rapid iteration.

## How It Works

1. The UI sends pathway inhibition updates (`0.0 -> 1.0`) to the backend.
2. The backend maps those controls to converted model parameters (`Ltot`, `Atot`, `FSK`, `IBMX`).
3. The converted ODE model advances by simulation steps.
4. The API returns a single state envelope with:
   - pathway activity values
   - ECG proxy trace points
   - pumping phase and pump metric

## Architecture

- `frontend/`: React + Vite app with interactive signaling diagram, ECG panel, and heart mechanics panel.
- `backend/`: FastAPI simulation API.
- `backend/converted_model/`: translated MATLAB equations and scenario helpers.
- `models/`: original MATLAB assets (`.m`, `.mat`, calibration inputs).

## API

- `GET /health`
- `GET /state`
- `POST /simulate/step`
- `POST /interventions`
- `POST /reset`

The frontend contract uses `SimulationEnvelope` and remains stable while model internals evolve.

## Local Run

### Backend

```bash
cd /Users/dev/Documents/hackathons/cerebral_valley/codex/backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

### Frontend

```bash
cd /Users/dev/Documents/hackathons/cerebral_valley/codex/frontend
npm install
npm run dev
```

Optional API base override:

```bash
VITE_API_BASE_URL=http://localhost:8000 npm run dev
```

## Demo Flow

1. Start from baseline simulation.
2. Click signaling nodes to apply inhibition levels.
3. Observe pathway state shifts in the network panel.
4. Observe electrical trace and pumping metric changes.
5. Reset and run a second intervention scenario.
