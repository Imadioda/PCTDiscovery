"""Text-only provider adapters. No embedded credentials, retries or redirects."""
from urllib.parse import urlparse, quote
import requests

PROVIDERS = {
    "Gemini": {"kind": "gemini", "base": "https://generativelanguage.googleapis.com/v1beta", "env": "GEMINI_API_KEY", "hint": "gemini-3.8-flash"},
    "OpenAI": {"kind": "chat", "base": "https://api.openai.com/v1", "env": "OPENAI_API_KEY", "hint": "Identifiant d'un modèle compatible Chat Completions", "token_parameter": "max_completion_tokens"},
    "Groq": {"kind": "chat", "base": "https://api.groq.com/openai/v1", "env": "GROQ_API_KEY", "hint": "openai/gpt-oss-120b", "token_parameter": "max_completion_tokens"},
    "Anthropic": {"kind": "anthropic", "base": "https://api.anthropic.com/v1", "env": "ANTHROPIC_API_KEY", "hint": "Identifiant Claude disponible dans ton compte"},
    "Mistral": {"kind": "chat", "base": "https://api.mistral.ai/v1", "env": "MISTRAL_API_KEY", "hint": "Identifiant disponible dans ton compte Mistral", "token_parameter": "max_tokens"},
    "OpenRouter": {"kind": "chat", "base": "https://openrouter.ai/api/v1", "env": "OPENROUTER_API_KEY", "hint": "fournisseur/identifiant-du-modele", "token_parameter": "max_tokens"},
    "xAI": {"kind": "chat", "base": "https://api.x.ai/v1", "env": "XAI_API_KEY", "hint": "Identifiant Grok disponible dans ton compte", "token_parameter": "max_tokens"},
    "Compatible OpenAI": {"kind": "chat", "base": "", "env": "", "hint": "Identifiant exact proposé par le service", "token_parameter": "max_tokens"},
}


class LLMError(RuntimeError):
    """Safe, user-facing error. Never includes credentials or raw response body."""
    def __init__(self, message, status=None):
        super().__init__(message)
        self.status = status


def validated_base_url(value):
    value = value.strip().rstrip("/")
    parsed = urlparse(value)
    local = parsed.hostname in {"localhost", "127.0.0.1", "::1"}
    if parsed.scheme != "https" and not (parsed.scheme == "http" and local):
        raise LLMError("Adresse HTTPS requise ; HTTP est accepté uniquement en local.")
    if not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise LLMError("Adresse de base invalide : pas d'identifiants, de paramètres ni de fragment dans l'URL.")
    if parsed.path.rstrip("/").endswith(("chat/completions", "messages", "generateContent")):
        raise LLMError("Saisis l'adresse de base (souvent terminée par /v1), pas l'endpoint complet.")
    return value


def status_error(status):
    messages = {
        400: "Requête refusée : vérifie le modèle, sa compatibilité et la taille du contexte.",
        401: "Clé API refusée par le fournisseur sélectionné.",
        402: "Crédit ou facturation insuffisants pour ce service.",
        403: "Accès refusé : vérifie les permissions de la clé et du modèle.",
        404: "Modèle ou endpoint introuvable : vérifie l'identifiant exact et le fournisseur.",
        413: "Contexte trop volumineux : réduis le nombre de gènes ou d'articles.",
        422: "Paramètres incompatibles avec ce modèle ou ce service.",
        429: "Quota ou limite de débit atteint. Consulte les limites du fournisseur avant de réessayer.",
        503: "Service temporairement indisponible ou surchargé. Réessaie plus tard ; ce statut ne prouve pas un dépassement de quota.",
        529: "Service surchargé. Réessaie plus tard.",
    }
    return LLMError(f"HTTP {status} — " + messages.get(status, "Le fournisseur n'a pas pu traiter la requête."), status)


def _post(url, headers, body):
    try:
        with requests.Session() as session:
            response = session.post(url, headers=headers, json=body,
                                    timeout=(15, 180), allow_redirects=False)
            with response:
                if not 200 <= response.status_code < 300:
                    raise status_error(response.status_code)
                try:
                    data = response.json()
                except ValueError:
                    raise LLMError("Le fournisseur a renvoyé une réponse non JSON.") from None
                if not isinstance(data, dict):
                    raise LLMError("Format de réponse non compatible.")
                if data.get("error"):
                    # Some compatible services wrap an error in HTTP 200.
                    code = data["error"].get("code") if isinstance(data["error"], dict) else None
                    if isinstance(code, int):
                        raise status_error(code)
                    raise LLMError("Le fournisseur a renvoyé une erreur dans sa réponse.")
                return data
    except requests.Timeout:
        raise LLMError("Le délai d'attente du fournisseur est dépassé. Aucun nouvel essai automatique.") from None
    except requests.RequestException:
        raise LLMError("Connexion au fournisseur impossible. Vérifie le réseau et l'adresse du service.") from None


def generate_text(provider, model, api_key, system_prompt, prompt,
                  base_url=None, max_output_tokens=8192,
                  token_parameter=None):
    if provider not in PROVIDERS:
        raise LLMError("Fournisseur non pris en charge.")
    config = PROVIDERS[provider]
    model = model.strip()
    if not model or len(model) > 256:
        raise LLMError("Indique l'identifiant exact du modèle (256 caractères maximum).")
    if not 512 <= max_output_tokens <= 16384:
        raise LLMError("La limite de sortie doit être comprise entre 512 et 16384 tokens.")
    custom = provider == "Compatible OpenAI"
    base = validated_base_url((base_url or "") if custom else config["base"])
    key = (api_key or "").strip()
    if not key and not custom:
        raise LLMError("Clé API manquante pour ce fournisseur.")
    headers = {"Content-Type": "application/json", "Accept-Encoding": "identity"}
    kind = config["kind"]
    if kind == "gemini":
        headers["x-goog-api-key"] = key
        url = base + "/models/" + quote(model.removeprefix("models/"), safe="") + ":generateContent"
        body = {"systemInstruction": {"parts": [{"text": system_prompt}]},
                "contents": [{"role": "user", "parts": [{"text": prompt}]}],
                "generationConfig": {"maxOutputTokens": max_output_tokens}}
    elif kind == "anthropic":
        headers.update({"x-api-key": key, "anthropic-version": "2023-06-01"})
        url = base + "/messages"
        body = {"model": model, "system": system_prompt, "max_tokens": max_output_tokens,
                "messages": [{"role": "user", "content": prompt}]}
    else:
        if key:
            headers["Authorization"] = f"Bearer {key}"
        url = base + "/chat/completions"
        parameter = token_parameter if custom else config["token_parameter"]
        parameter = parameter or "max_tokens"
        if parameter not in {"max_tokens", "max_completion_tokens"}:
            raise LLMError("Paramètre de limite de sortie non pris en charge.")
        body = {"model": model, "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": prompt}],
            parameter: max_output_tokens, "stream": False}
        if provider == "OpenAI":
            body["store"] = False
    data = _post(url, headers, body)
    try:
        if kind == "gemini":
            candidates = data.get("candidates") or []
            first = candidates[0] if candidates else {}
            text = "".join(part.get("text", "") for part in first.get("content", {}).get("parts", []) if not part.get("thought"))
            reason = first.get("finishReason", "NO_CANDIDATE")
            usage = data.get("usageMetadata", {})
            returned = data.get("modelVersion", model)
            complete = reason == "STOP"
        elif kind == "anthropic":
            text = "".join(part.get("text", "") for part in data.get("content", []) if part.get("type") == "text")
            reason = data.get("stop_reason", "UNKNOWN")
            usage = data.get("usage", {})
            returned = data.get("model", model)
            complete = reason == "end_turn"
        else:
            choices = data.get("choices") or []
            first = choices[0] if choices else {}
            content = first.get("message", {}).get("content")
            if isinstance(content, list):
                text = "".join(part.get("text", "") for part in content if part.get("type") == "text")
            else:
                text = content or ""
            reason = first.get("finish_reason", "NO_CHOICE")
            usage = data.get("usage", {})
            returned = data.get("model", model)
            complete = reason == "stop"
        if not isinstance(text, str):
            raise TypeError()
    except (TypeError, AttributeError, IndexError):
        raise LLMError("La structure de la réponse n'est pas compatible avec cet adaptateur.") from None
    return {"text": text, "finish_reason": reason, "complete": complete,
            "usage": usage, "model": returned}
