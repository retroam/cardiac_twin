function traceToPolylinePoints(trace, width = 560, height = 210, padding = 14) {
  if (!trace || trace.length < 2) {
    return "";
  }

  return trace
    .map((value, index) => {
      const x = padding + (index / (trace.length - 1)) * (width - padding * 2);
      const y = height / 2 - value * ((height - padding * 2) / 2);
      return `${x.toFixed(2)},${y.toFixed(2)}`;
    })
    .join(" ");
}

export default function ECGPanel({ tracePoints, heartRate }) {
  const points = traceToPolylinePoints(tracePoints);

  return (
    <section className="panel ecg-panel">
      <div className="panel-header">
        <h2>Electrical Trace (ECG Proxy)</h2>
        <span>{heartRate.toFixed(1)} BPM</span>
      </div>

      <svg viewBox="0 0 560 210" role="img" aria-label="ECG waveform">
        <defs>
          <pattern id="ecg-grid" width="20" height="20" patternUnits="userSpaceOnUse">
            <path d="M 20 0 L 0 0 0 20" fill="none" stroke="rgba(116,135,153,0.24)" strokeWidth="1" />
          </pattern>
        </defs>
        <rect x="0" y="0" width="560" height="210" fill="url(#ecg-grid)" />
        <line x1="0" y1="105" x2="560" y2="105" stroke="rgba(64,82,100,0.25)" strokeWidth="1.5" />
        <polyline points={points} fill="none" stroke="var(--accent)" strokeWidth="3" strokeLinecap="round" />
      </svg>
    </section>
  );
}
