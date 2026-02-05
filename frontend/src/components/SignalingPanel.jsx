import { useMemo } from "react";

const PROTEIN_NODES = [
  {
    id: "b1ar",
    label: "b1AR",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 470,
    y: 76,
    rx: 44,
    ry: 30
  },
  {
    id: "bark",
    label: "bARK",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 340,
    y: 106,
    rx: 36,
    ry: 23
  },
  {
    id: "gsa_l",
    label: "Gsa",
    pathwayId: "beta_adrenergic",
    shape: "triangle",
    x: 526,
    y: 106,
    w: 66,
    h: 56,
    dy: 4
  },
  {
    id: "bg_l",
    label: "bg",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 495,
    y: 144,
    rx: 31,
    ry: 19
  },
  {
    id: "gsa_r",
    label: "Gsa",
    pathwayId: "beta_adrenergic",
    shape: "triangle",
    x: 604,
    y: 108,
    w: 66,
    h: 56,
    dy: 4
  },
  {
    id: "bg_r",
    label: "bg",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 566,
    y: 146,
    rx: 31,
    ry: 19
  },
  {
    id: "ac",
    label: "AC",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 690,
    y: 88,
    rx: 42,
    ry: 28
  },
  {
    id: "pka_top",
    label: "PKA",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 579,
    y: 205,
    rx: 36,
    ry: 24
  },
  {
    id: "pde3",
    label: "PDE3",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 736,
    y: 208,
    rx: 42,
    ry: 24
  },
  {
    id: "pde4",
    label: "PDE4",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 816,
    y: 188,
    rx: 42,
    ry: 24
  },
  {
    id: "ibmx",
    label: "IBMX",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 806,
    y: 271,
    rx: 45,
    ry: 27
  },
  {
    id: "pki",
    label: "PKI",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 449,
    y: 270,
    rx: 34,
    ry: 22
  },
  {
    id: "pka_mid",
    label: "PKA",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 533,
    y: 270,
    rx: 36,
    ry: 24
  },
  {
    id: "r_site",
    label: "r",
    shape: "ellipse",
    x: 620,
    y: 267,
    rx: 16,
    ry: 14,
    dy: 5,
    interactive: false
  },
  {
    id: "i1",
    label: "I1",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 606,
    y: 314,
    rx: 32,
    ry: 22
  },
  {
    id: "pp2a",
    label: "PP2A",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 802,
    y: 314,
    rx: 44,
    ry: 24
  },
  {
    id: "pp1",
    label: "PP1",
    pathwayId: "contractility",
    shape: "ellipse",
    x: 540,
    y: 375,
    rx: 36,
    ry: 24
  },
  {
    id: "plb",
    label: "PLB",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 465,
    y: 317,
    rx: 34,
    ry: 22
  },
  {
    id: "serca",
    label: "SERCA",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 377,
    y: 330,
    rx: 47,
    ry: 24
  },
  {
    id: "ryr",
    label: "RyR",
    pathwayId: "calcium_handling",
    shape: "rect",
    x: 327,
    y: 293,
    w: 80,
    h: 52,
    r: 4
  },
  {
    id: "cam",
    label: "CaM",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 360,
    y: 426,
    rx: 36,
    ry: 24
  },
  {
    id: "tnc",
    label: "TnC",
    pathwayId: "contractility",
    shape: "ellipse",
    x: 646,
    y: 423,
    rx: 32,
    ry: 22
  },
  {
    id: "tni",
    label: "TnI",
    pathwayId: "contractility",
    shape: "ellipse",
    x: 775,
    y: 423,
    rx: 32,
    ry: 22
  }
];

const ACTIVATION_EDGES = [
  "M470 20 L470 44",
  "M374 104 C405 100 431 91 448 83",
  "M503 88 C512 91 519 96 523 101",
  "M482 102 C488 117 492 127 495 138",
  "M559 108 L571 108",
  "M637 104 C651 98 668 93 681 90",
  "M710 65 L724 76",
  "M688 116 C680 136 670 148 657 160",
  "M648 170 C632 184 610 194 591 202",
  "M575 228 C568 242 558 255 547 264",
  "M589 226 C598 248 604 274 606 292",
  "M563 282 C576 293 587 302 596 308",
  "M506 281 C494 294 483 304 474 312",
  "M560 284 C624 330 693 369 753 409",
  "M394 424 C480 413 565 412 613 420",
  "M416 323 C392 314 366 305 344 298",
  "M332 319 C338 356 347 386 354 404",
  "M612 223 C621 238 628 247 635 255"
];

const INHIBITION_EDGES = [
  {
    d: "M723 201 C704 188 682 176 664 168",
    bar: { x1: 657, y1: 160, x2: 670, y2: 174 }
  },
  {
    d: "M812 185 C757 164 702 152 664 151",
    bar: { x1: 656, y1: 143, x2: 672, y2: 158 }
  },
  {
    d: "M805 244 C810 226 815 211 817 196",
    bar: { x1: 809, y1: 194, x2: 825, y2: 198 }
  },
  {
    d: "M789 258 C770 241 754 228 744 219",
    bar: { x1: 738, y1: 213, x2: 752, y2: 222 }
  },
  {
    d: "M482 270 L507 270",
    bar: { x1: 509, y1: 262, x2: 509, y2: 278 }
  },
  {
    d: "M588 336 C575 352 560 365 548 372",
    bar: { x1: 541, y1: 366, x2: 554, y2: 378 }
  },
  {
    d: "M757 314 L640 314",
    bar: { x1: 637, y1: 306, x2: 637, y2: 322 }
  },
  {
    d: "M509 360 C492 343 480 330 470 322",
    bar: { x1: 465, y1: 315, x2: 476, y2: 326 }
  },
  {
    d: "M431 321 C417 322 405 324 391 326",
    bar: { x1: 388, y1: 318, x2: 388, y2: 334 }
  },
  {
    d: "M570 385 C648 400 705 412 744 420",
    bar: { x1: 746, y1: 410, x2: 746, y2: 430 }
  },
  {
    d: "M678 423 L742 423",
    bar: { x1: 744, y1: 414, x2: 744, y2: 432 }
  }
];

const RED_ANNOTATIONS = [
  { text: "ADRBK1", x: 245, y: 137 },
  { text: "ADRBK2", x: 245, y: 161 },
  { text: "GNAS", x: 541, y: 140 },
  { text: "PRKACA", x: 420, y: 171 },
  { text: "PRKACB", x: 420, y: 192 },
  { text: "PRKAR1A", x: 420, y: 213 },
  { text: "PRKAR1B", x: 420, y: 234 },
  { text: "PDE3A", x: 658, y: 230 },
  { text: "PDE3B", x: 658, y: 252 },
  { text: "PDE4A", x: 860, y: 196 },
  { text: "PDE4B", x: 860, y: 218 },
  { text: "PDE4C", x: 860, y: 240 },
  { text: "PDE4D", x: 860, y: 262 },
  { text: "PP1R1B", x: 629, y: 291 },
  { text: "PP1R1C", x: 629, y: 313 },
  { text: "PP1R1A", x: 629, y: 350 },
  { text: "PPP1CA", x: 427, y: 394 },
  { text: "PPP1CB", x: 427, y: 416 },
  { text: "PPP1CC", x: 427, y: 438 },
  { text: "PLN", x: 498, y: 336 },
  { text: "ATP2A2", x: 378, y: 357 },
  { text: "RyR2", x: 310, y: 282 },
  { text: "CALM1", x: 412, y: 427 },
  { text: "CALM2", x: 412, y: 449 },
  { text: "CALM3", x: 412, y: 471 },
  { text: "TNNI3", x: 828, y: 425 }
];

function toPercent(value) {
  return `${Math.round(value * 100)}%`;
}

function nextInhibitionLevel(current) {
  if (current < 0.2) {
    return 0.5;
  }
  if (current < 0.7) {
    return 0.85;
  }
  return 0;
}

function inhibitionClass(inhibition) {
  if (inhibition >= 0.7) {
    return "blocked";
  }
  if (inhibition >= 0.2) {
    return "damped";
  }
  return "normal";
}

function renderNodeShape(node) {
  if (node.shape === "triangle") {
    const halfW = node.w / 2;
    const halfH = node.h / 2;
    const points = `${node.x - halfW},${node.y - halfH} ${node.x + halfW},${node.y - halfH} ${node.x},${node.y + halfH}`;
    return <polygon points={points} className="protein-shape" />;
  }

  if (node.shape === "rect") {
    return (
      <rect
        x={node.x - node.w / 2}
        y={node.y - node.h / 2}
        width={node.w}
        height={node.h}
        rx={node.r || 0}
        className="protein-shape"
      />
    );
  }

  return <ellipse cx={node.x} cy={node.y} rx={node.rx} ry={node.ry} className="protein-shape" />;
}

export default function SignalingPanel({
  pathways,
  onInhibitionChange,
  onToggleRun,
  onReset,
  isRunning,
  busy
}) {
  const pathwayById = useMemo(
    () =>
      pathways.reduce((accumulator, pathway) => {
        accumulator[pathway.id] = pathway;
        return accumulator;
      }, {}),
    [pathways]
  );

  return (
    <section className="panel signaling-panel">
      <div className="panel-header">
        <h2>Signaling Cartoon (Reference-Matched)</h2>
        <div className="controls-actions">
          <button type="button" className="ghost-button" onClick={onToggleRun}>
            {isRunning ? "Pause" : "Resume"}
          </button>
          <button type="button" className="ghost-button" onClick={onReset} disabled={busy}>
            Reset
          </button>
        </div>
      </div>

      <p className="panel-copy">Click proteins directly in the pathway map to adjust inhibition levels.</p>

      <div className="network-canvas-wrap">
        <svg className="network-canvas" viewBox="0 0 1060 580" role="img" aria-label="Cardiac signaling proteins and arrows">
          <defs>
            <marker id="arrowhead" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto">
              <path d="M0,0 L8,4 L0,8 z" fill="#111" />
            </marker>
          </defs>

          <rect x="20" y="18" width="1020" height="552" rx="60" className="cell-backdrop" />
          <path
            className="ttubule-path"
            d="M130,48 C88,50 78,96 84,142 C90,188 124,202 148,205 C174,208 189,231 177,256 C166,279 138,305 146,347 C154,391 189,472 255,470"
          />

          <path
            className="sr-outline"
            d="M312 248 C280 258 272 301 302 341 C325 372 372 365 404 338 C428 318 432 285 408 268 C379 249 336 254 317 281"
          />

          <path className="myo-outline" d="M602 402 L830 402 L874 438 L830 474 L602 474 L558 438 Z" />
          <path className="myo-rail" d="M620 417 L814 417" />
          <path className="myo-rail" d="M613 432 L806 432" />
          <path className="myo-rail" d="M620 447 L814 447" />

          <text x="208" y="487" className="region-label">T-tubule</text>
          <text x="626" y="508" className="region-label">Myofilaments</text>

          <text x="418" y="23" className="annotation-text">NE, ISO</text>
          <text x="721" y="52" className="annotation-text">Fsk</text>
          <text x="756" y="100" className="annotation-text">ATP5&apos;-AMP</text>
          <text x="648" y="168" className="annotation-text">cAMP</text>
          <text x="329" y="350" className="annotation-text sr-label">SR</text>

          {ACTIVATION_EDGES.map((pathD, index) => (
            <path key={`act-${index}`} d={pathD} className="edge-activation" markerEnd="url(#arrowhead)" />
          ))}

          {INHIBITION_EDGES.map((edge, index) => (
            <g key={`inh-${index}`}>
              <path d={edge.d} className="edge-inhibition" />
              <line x1={edge.bar.x1} y1={edge.bar.y1} x2={edge.bar.x2} y2={edge.bar.y2} className="edge-tbar" />
            </g>
          ))}

          {RED_ANNOTATIONS.map((item) => (
            <text key={item.text} x={item.x} y={item.y} className="gene-label">
              {item.text}
            </text>
          ))}

          {PROTEIN_NODES.map((node) => {
            const pathway = node.pathwayId ? pathwayById[node.pathwayId] : null;
            const inhibition = pathway?.inhibition ?? 0;
            const nodeClass = inhibitionClass(inhibition);
            const isInteractive = node.interactive !== false && Boolean(node.pathwayId);

            return (
              <g
                key={node.id}
                className={`protein-node ${isInteractive ? "interactive" : "static"} ${nodeClass}`}
                role={isInteractive ? "button" : undefined}
                tabIndex={isInteractive ? 0 : -1}
                aria-label={isInteractive ? `${node.label} toggle inhibition` : undefined}
                onClick={
                  isInteractive
                    ? () => onInhibitionChange({ [node.pathwayId]: nextInhibitionLevel(inhibition) })
                    : undefined
                }
                onKeyDown={
                  isInteractive
                    ? (event) => {
                        if (event.key === "Enter" || event.key === " ") {
                          event.preventDefault();
                          onInhibitionChange({ [node.pathwayId]: nextInhibitionLevel(inhibition) });
                        }
                      }
                    : undefined
                }
              >
                {renderNodeShape(node)}
                <text x={node.x} y={node.y + (node.dy ?? 6)} textAnchor="middle" className="protein-label">
                  {node.label}
                </text>
              </g>
            );
          })}
        </svg>
      </div>

      <div className="legend-row">
        <div className="legend-item">
          <span className="legend-swatch normal" />
          Normal
        </div>
        <div className="legend-item">
          <span className="legend-swatch damped" />
          Partial inhibition
        </div>
        <div className="legend-item">
          <span className="legend-swatch blocked" />
          Strong inhibition
        </div>
      </div>

      <div className="pathway-readout">
        {pathways.map((pathway) => (
          <article className="pathway-chip" key={pathway.id}>
            <h3>{pathway.label}</h3>
            <p>{toPercent(pathway.inhibition)} inhibition</p>
            <p>{toPercent(pathway.activity)} activity</p>
          </article>
        ))}
      </div>
    </section>
  );
}
