import { useEffect, useState } from "react";

const PRESETS = [
  {
    id: "beta_blockade",
    label: "Beta Blockade",
    pathwayIds: ["beta_adrenergic"],
    strength: 0.72
  },
  {
    id: "calcium_suppression",
    label: "Calcium Suppression",
    pathwayIds: ["calcium_handling"],
    strength: 0.64
  },
  {
    id: "combined_stress",
    label: "Combined Stress",
    pathwayIds: ["beta_adrenergic", "ion_channel_current"],
    strength: 0.48
  }
];

function toPercent(value) {
  return `${Math.round(value * 100)}%`;
}

export default function PathwayControls({
  pathways,
  onApplyIntervention,
  onApplyPreset,
  onToggleRun,
  onReset,
  isRunning,
  busy
}) {
  const [draftInhibition, setDraftInhibition] = useState({});

  useEffect(() => {
    const next = {};
    pathways.forEach((pathway) => {
      next[pathway.id] = pathway.inhibition;
    });
    setDraftInhibition(next);
  }, [pathways]);

  return (
    <section className="panel controls-panel">
      <div className="panel-header">
        <h2>Pathway Interventions</h2>
        <div className="controls-actions">
          <button type="button" className="ghost-button" onClick={onToggleRun}>
            {isRunning ? "Pause" : "Resume"}
          </button>
          <button type="button" className="ghost-button" onClick={onReset} disabled={busy}>
            Reset
          </button>
        </div>
      </div>

      <p className="panel-copy">
        Each slider controls inhibition strength for one pathway (0.00 to 1.00).
      </p>

      <div className="pathway-grid">
        {pathways.map((pathway) => {
          const draftValue = draftInhibition[pathway.id] ?? pathway.inhibition;
          return (
            <div className="pathway-row" key={pathway.id}>
              <div className="pathway-title">
                <h3>{pathway.label}</h3>
                <span>{toPercent(pathway.activity)} active</span>
              </div>

              <input
                aria-label={`${pathway.label} inhibition`}
                type="range"
                min="0"
                max="1"
                step="0.01"
                value={draftValue}
                onChange={(event) => {
                  const value = Number(event.target.value);
                  setDraftInhibition((prev) => ({ ...prev, [pathway.id]: value }));
                }}
              />

              <div className="pathway-actions">
                <span>Inhibit {toPercent(draftValue)}</span>
                <button
                  type="button"
                  onClick={() => onApplyIntervention([pathway.id], draftValue)}
                  disabled={busy}
                >
                  Apply
                </button>
              </div>
            </div>
          );
        })}
      </div>

      <div className="preset-strip">
        <span>Demo presets:</span>
        {PRESETS.map((preset) => (
          <button
            type="button"
            key={preset.id}
            className="preset-button"
            onClick={() => onApplyPreset(preset.pathwayIds, preset.strength)}
            disabled={busy}
          >
            {preset.label}
          </button>
        ))}
      </div>
    </section>
  );
}
