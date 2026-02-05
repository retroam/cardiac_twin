function toPercent(value) {
  return `${Math.round(value * 100)}%`;
}

export default function HeartPanel({ mechanical, heartRate }) {
  const oscillation = 0.5 + 0.5 * Math.sin(mechanical.phase_radians);
  const beatScale = 0.82 + mechanical.pump_metric * (0.22 + 0.1 * oscillation);

  return (
    <section className="panel heart-panel">
      <div className="panel-header">
        <h2>Pumping Mechanics</h2>
      </div>

      <div className="heart-stage">
        <div
          className={`heart-core ${mechanical.phase}`}
          style={{ transform: `translate(-50%, -50%) scale(${beatScale.toFixed(3)})` }}
          aria-hidden="true"
        />
      </div>

      <div className="metric-row">
        <span>Pump Metric (stroke-volume proxy)</span>
        <strong>{toPercent(mechanical.pump_metric)}</strong>
      </div>
      <div className="meter-track">
        <div className="meter-fill" style={{ width: `${Math.max(mechanical.pump_metric * 100, 2)}%` }} />
      </div>

      <div className="stat-grid">
        <div>
          <span>Cycle Phase</span>
          <strong>{mechanical.phase}</strong>
        </div>
        <div>
          <span>Heart Rate</span>
          <strong>{heartRate.toFixed(1)} BPM</strong>
        </div>
      </div>
    </section>
  );
}
