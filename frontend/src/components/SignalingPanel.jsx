import { useMemo } from "react";

const PROTEIN_NODES = [
  {
    id: "b1ar",
    label: "b1AR",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 472,
    y: 78,
    rx: 44,
    ry: 30
  },
  {
    id: "bark",
    label: "bARK",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 342,
    y: 110,
    rx: 36,
    ry: 23
  },
  {
    id: "gsa_l",
    label: "Gsa",
    pathwayId: "beta_adrenergic",
    shape: "triangle",
    x: 528,
    y: 108,
    w: 66,
    h: 56,
    dy: 4
  },
  {
    id: "bg_l",
    label: "bg",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 498,
    y: 146,
    rx: 31,
    ry: 19
  },
  {
    id: "gsa_r",
    label: "Gsa",
    pathwayId: "beta_adrenergic",
    shape: "triangle",
    x: 611,
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
    x: 572,
    y: 146,
    rx: 31,
    ry: 19
  },
  {
    id: "ac",
    label: "AC",
    pathwayId: "beta_adrenergic",
    shape: "ellipse",
    x: 696,
    y: 90,
    rx: 42,
    ry: 28
  },
  {
    id: "pka_top",
    label: "PKA",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 578,
    y: 214,
    rx: 36,
    ry: 24
  },
  {
    id: "pde3",
    label: "PDE3",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 736,
    y: 216,
    rx: 42,
    ry: 24
  },
  {
    id: "pde4",
    label: "PDE4",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 815,
    y: 190,
    rx: 42,
    ry: 24
  },
  {
    id: "ibmx",
    label: "IBMX",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 804,
    y: 276,
    rx: 45,
    ry: 27
  },
  {
    id: "pki",
    label: "PKI",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 466,
    y: 277,
    rx: 34,
    ry: 22
  },
  {
    id: "pka_mid",
    label: "PKA",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 536,
    y: 277,
    rx: 36,
    ry: 24
  },
  {
    id: "r_site",
    label: "r",
    shape: "ellipse",
    x: 620,
    y: 273,
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
    y: 330,
    rx: 32,
    ry: 22
  },
  {
    id: "pp2a",
    label: "PP2A",
    pathwayId: "ion_channel_current",
    shape: "ellipse",
    x: 802,
    y: 330,
    rx: 44,
    ry: 24
  },
  {
    id: "pp1",
    label: "PP1",
    pathwayId: "contractility",
    shape: "ellipse",
    x: 540,
    y: 380,
    rx: 36,
    ry: 24
  },
  {
    id: "plb",
    label: "PLB",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 468,
    y: 323,
    rx: 34,
    ry: 22
  },
  {
    id: "serca",
    label: "SERCA",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 382,
    y: 334,
    rx: 47,
    ry: 24
  },
  {
    id: "ryr",
    label: "RyR",
    pathwayId: "calcium_handling",
    shape: "rect",
    x: 322,
    y: 296,
    w: 80,
    h: 52,
    r: 4
  },
  {
    id: "cam",
    label: "CaM",
    pathwayId: "calcium_handling",
    shape: "ellipse",
    x: 356,
    y: 436,
    rx: 36,
    ry: 24
  },
  {
    id: "tnc",
    label: "TnC",
    pathwayId: "contractility",
    shape: "ellipse",
    x: 648,
    y: 431,
    rx: 32,
    ry: 22
  },
  {
    id: "tni",
    label: "TnI",
    pathwayId: "contractility",
    shape: "ellipse",
    x: 774,
    y: 431,
    rx: 32,
    ry: 22
  }
];

const ACTIVATION_EDGES = [
  "M472 22 L472 47",
  "M376 110 C408 103 436 95 452 86",
  "M505 92 C514 95 521 101 525 106",
  "M486 104 C490 119 494 130 498 140",
  "M562 108 L580 108",
  "M644 106 C658 98 678 94 689 92",
  "M725 58 L739 78",
  "M695 120 C684 144 672 160 654 172",
  "M650 181 C633 196 607 208 589 212",
  "M674 176 C695 189 712 202 726 212",
  "M677 170 C717 168 752 173 785 184",
  "M575 237 C568 252 556 266 544 273",
  "M566 288 C578 302 591 312 598 322",
  "M508 286 C495 301 482 311 474 320",
  "M638 340 C678 365 718 392 754 424",
  "M402 433 C486 425 571 424 618 428",
  "M418 329 C393 320 366 310 339 300",
  "M332 321 C337 364 348 401 352 412",
  "M610 232 C620 248 628 259 637 266"
];

const INHIBITION_EDGES = [
  {
    d: "M804 248 C808 229 814 211 815 198",
    bar: { x1: 807, y1: 196, x2: 823, y2: 200 }
  },
  {
    d: "M786 264 C769 246 753 232 744 223",
    bar: { x1: 738, y1: 217, x2: 752, y2: 226 }
  },
  {
    d: "M500 277 L514 277",
    bar: { x1: 517, y1: 269, x2: 517, y2: 285 }
  },
  {
    d: "M594 350 C580 365 565 375 550 379",
    bar: { x1: 543, y1: 371, x2: 556, y2: 385 }
  },
  {
    d: "M758 330 L642 330",
    bar: { x1: 639, y1: 322, x2: 639, y2: 338 }
  },
  {
    d: "M510 367 C495 349 482 336 474 327",
    bar: { x1: 468, y1: 319, x2: 479, y2: 330 }
  },
  {
    d: "M434 327 C418 328 406 329 394 331",
    bar: { x1: 391, y1: 323, x2: 391, y2: 339 }
  },
  {
    d: "M570 388 C650 402 705 416 744 428",
    bar: { x1: 746, y1: 420, x2: 746, y2: 436 }
  },
  {
    d: "M801 355 C798 382 792 406 783 423",
    bar: { x1: 775, y1: 421, x2: 790, y2: 429 }
  },
  {
    d: "M681 431 L742 431",
    bar: { x1: 744, y1: 423, x2: 744, y2: 439 }
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
  { text: "PP1R1A", x: 628, y: 370 },
  { text: "PPP1CA", x: 428, y: 401 },
  { text: "PPP1CB", x: 428, y: 423 },
  { text: "PPP1CC", x: 428, y: 445 },
  { text: "PLN", x: 498, y: 346 },
  { text: "ATP2A2", x: 380, y: 365 },
  { text: "RyR2", x: 310, y: 282 },
  { text: "CALM1", x: 410, y: 437 },
  { text: "CALM2", x: 410, y: 459 },
  { text: "CALM3", x: 410, y: 481 },
  { text: "TNNI3", x: 834, y: 430 }
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

          <text x="424" y="26" className="annotation-text">NE, ISO</text>
          <text x="716" y="54" className="annotation-text">Fsk</text>
          <text x="754" y="102" className="annotation-text">ATP5&apos;-AMP</text>
          <text x="648" y="170" className="annotation-text">cAMP</text>
          <text x="328" y="356" className="annotation-text sr-label">SR</text>

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
