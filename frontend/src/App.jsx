import { useCallback, useEffect, useMemo, useRef, useState } from "react";

import { API_BASE_URL, applyInterventions, getState, resetSimulation, stepSimulation } from "./api/client";
import ECGPanel from "./components/ECGPanel";
import HeartPanel from "./components/HeartPanel";
import SignalingPanel from "./components/SignalingPanel";

async function withFallback(task, fallbackMessage) {
  try {
    return await task();
  } catch (error) {
    if (error instanceof Error) {
      throw error;
    }
    throw new Error(fallbackMessage);
  }
}

export default function App() {
  const [envelope, setEnvelope] = useState(null);
  const [error, setError] = useState("");
  const [busy, setBusy] = useState(false);
  const [isRunning, setIsRunning] = useState(true);
  const [statusText, setStatusText] = useState("Connecting to simulation backend...");

  const inFlightRef = useRef(false);

  const runRequest = useCallback(async (requestFactory, nextStatus = "Synced") => {
    if (inFlightRef.current) {
      return;
    }

    inFlightRef.current = true;
    setBusy(true);

    try {
      const payload = await requestFactory();
      setEnvelope(payload);
      setError("");
      setStatusText(nextStatus);
    } catch (requestError) {
      setError(requestError instanceof Error ? requestError.message : "Unknown request error");
      setStatusText("Backend connection failed");
    } finally {
      inFlightRef.current = false;
      setBusy(false);
    }
  }, []);

  useEffect(() => {
    runRequest(
      () => withFallback(() => getState(), "Failed to fetch simulation state."),
      "Ready"
    );
  }, [runRequest]);

  useEffect(() => {
    if (!isRunning) {
      return undefined;
    }

    const timer = setInterval(() => {
      runRequest(
        () => withFallback(() => stepSimulation({ steps: 1, dt_ms: 70 }), "Step failed."),
        "Live stepping"
      );
    }, 360);

    return () => clearInterval(timer);
  }, [isRunning, runRequest]);

  const state = envelope?.state;

  const pathways = useMemo(() => {
    if (!state) {
      return [];
    }
    return Object.values(state.pathways);
  }, [state]);

  const handleSetPathwayInhibitions = useCallback(
    (inhibitionByPathway) => {
      const entries = Object.entries(inhibitionByPathway);
      if (!entries.length) {
        return;
      }

      runRequest(
        async () => {
          let latestPayload = null;
          for (const [pathwayId, inhibitionStrength] of entries) {
            latestPayload = await applyInterventions({
              pathway_ids: [pathwayId],
              inhibition_strength: Number(inhibitionStrength.toFixed(2))
            });
          }

          if (!latestPayload) {
            return withFallback(() => getState(), "Failed to refresh simulation state.");
          }

          return latestPayload;
        },
        "Updated signaling network"
      );
    },
    [runRequest]
  );

  const handleReset = useCallback(() => {
    runRequest(() => withFallback(() => resetSimulation(42), "Reset failed."), "Reset to baseline");
  }, [runRequest]);

  if (!state) {
    return (
      <main className="app-shell">
        <section className="hero-card">
          <h1>CardioTwin</h1>
          <p>{statusText}</p>
          {error ? <p className="error-copy">{error}</p> : null}
        </section>
      </main>
    );
  }

  return (
    <main className="app-shell">
      <header className="hero-card animate-rise">
        <div>
          <p className="eyebrow">Transcriptome + Cardiac Systems Sandbox</p>
          <h1>CardioTwin Interactive Simulation</h1>
          <p className="subtitle">
            Click the signaling cartoon to perturb pathways and watch ECG and pumping changes propagate in real time.
          </p>
        </div>

        <div className="hero-stats">
          <div>
            <span>Backend</span>
            <strong>{envelope.meta.placeholder_backend ? "Placeholder" : "Live"}</strong>
          </div>
          <div>
            <span>Step</span>
            <strong>{state.step_index}</strong>
          </div>
          <div>
            <span>Time</span>
            <strong>{(state.time_ms / 1000).toFixed(2)}s</strong>
          </div>
          <div>
            <span>Status</span>
            <strong>{busy ? "Updating" : statusText}</strong>
          </div>
        </div>
      </header>

      {error ? <p className="error-copy">{error}</p> : null}

      <div className="layout-grid">
        <SignalingPanel
          pathways={pathways}
          onInhibitionChange={handleSetPathwayInhibitions}
          onToggleRun={() => setIsRunning((previous) => !previous)}
          onReset={handleReset}
          isRunning={isRunning}
          busy={busy}
        />

        <ECGPanel
          tracePoints={state.electrical.trace_points}
          heartRate={state.electrical.heart_rate_bpm}
        />

        <HeartPanel
          mechanical={state.mechanical}
          heartRate={state.electrical.heart_rate_bpm}
        />
      </div>

      <footer className="footer-note">
        <span>API Base:</span> {API_BASE_URL}
      </footer>
    </main>
  );
}
