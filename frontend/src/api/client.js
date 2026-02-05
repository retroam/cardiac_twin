const API_BASE_URL = import.meta.env.VITE_API_BASE_URL || "http://localhost:8000";

async function request(path, options = {}) {
  const response = await fetch(`${API_BASE_URL}${path}`, {
    headers: { "Content-Type": "application/json", ...(options.headers || {}) },
    ...options
  });

  if (!response.ok) {
    const text = await response.text();
    throw new Error(text || `Request failed (${response.status})`);
  }

  return response.json();
}

export async function getState() {
  return request("/state");
}

export async function stepSimulation(payload = { steps: 1, dt_ms: 60 }) {
  return request("/simulate/step", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export async function applyInterventions(payload) {
  return request("/interventions", {
    method: "POST",
    body: JSON.stringify(payload)
  });
}

export async function resetSimulation(seed = 42) {
  return request("/reset", {
    method: "POST",
    body: JSON.stringify({ seed })
  });
}

export { API_BASE_URL };
