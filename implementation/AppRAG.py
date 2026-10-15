"""Local interface. Run: python -m streamlit run AppRAG.py"""
from pathlib import Path
from datetime import datetime, timezone
import json
import re
from urllib.parse import urlparse
from hashlib import sha256

import pandas as pd
import streamlit as st

from AskRAG import answer_question, TOPICS, MODEL
from LLMProviders import PROVIDERS, LLMError, validated_base_url
from BuildRAGContext import RAGContextBuilder

BASE_DIR = Path(__file__).resolve().parent
CANCER_NAMES = {
    "BRCA_breast": "Sein — carcinome invasif",
    "COAD_colon": "Côlon — adénocarcinome",
    "KIRC_kidney": "Rein — carcinome à cellules claires",
    "KIRP_kidney": "Rein — carcinome papillaire",
    "LIHC_liver": "Foie — carcinome hépatocellulaire",
    "LUAD_lung": "Poumon — adénocarcinome",
    "LUSC_lung": "Poumon — carcinome épidermoïde",
    "PRAD_prostate": "Prostate — adénocarcinome",
    "STAD_stomach": "Estomac — adénocarcinome",
    "THCA_thyroid": "Thyroïde — carcinome",
    "UCEC_endometrium": "Endomètre — carcinome",
}
DIRECTIONS = {
    "Les deux directions": "both",
    "Moins exprimés dans les tumeurs": "down",
    "Plus exprimés dans les tumeurs": "up",
}
STATUS_NAMES = {
    "candidate_documents_retrieved": "Résumés disponibles — association à vérifier",
    "gene_mapping_requires_review": "Identifiant à vérifier",
    "not_searched_in_current_collection": "Pas encore recherché dans la collection",
    "collection_request_failed": "Échec de la collecte documentaire",
    "no_articles_found_by_collection_query": "Aucun article trouvé par notre requête",
    "no_documents_retrieved_after_filtering": "Aucun résumé retenu après filtrage",
}


def suggest_question(kind):
    prompts = {
        "explain": "Explique simplement les résultats des gènes sélectionnés et ce que les sources rapportent sur leur rôle.",
        "targets": "Quelles pistes thérapeutiques sont évoquées pour les gènes sélectionnés ? Distingue hypothèses, résultats précliniques et preuves cliniques, sans faire de recommandation de traitement.",
        "limits": "Quelles conclusions peut-on tirer de ces résultats, et quelles conclusions ne sont pas justifiées par les données et les sources ?",
    }
    st.session_state.question = prompts[kind]
    if kind == "targets":
        st.session_state.topic = "2"


def render_answer(entry, index):
    context, result = entry["context"], entry["result"]
    with st.chat_message("user"):
        st.write(entry["question"])
        st.caption(CANCER_NAMES.get(context["comparison_id"], context["comparison_id"]))
    with st.chat_message("assistant"):
        audit = result["audit"]
        st.caption(f"{audit.get('provider', '')} · {audit.get('model_returned', '')}")
        if audit["warnings"]:
            st.warning("Cette réponse nécessite une vérification : " + "; ".join(audit["warnings"]))
        st.markdown(audit["answer"], unsafe_allow_html=False)
        with st.expander("Sources documentaires et résumés"):
            ids = set(re.findall(r"\bMED:\d+\b", audit["answer"]))
            if not context["sources"]:
                st.info("Aucun résumé fourni pour cette sélection. Cela ne prouve pas une absence de publications.")
            for source in context["sources"]:
                title = re.sub(r"<[^>]*>", "", source["title"])
                st.write(f"{source['source_id']} · {title}")
                st.caption(f"Année : {source['year']} · " + ("Cité dans la réponse" if source['source_id'] in ids else "Fourni au modèle, non cité"))
                url = source.get("url", "")
                if urlparse(url).scheme == "https" and urlparse(url).hostname == "europepmc.org":
                    st.link_button("Consulter sur Europe PMC", url)
                st.write(source["abstract"])
                st.divider()
        with st.expander("Consommation de cet appel"):
            st.caption("Compteurs renvoyés par le fournisseur, pas le cumul du compte. Le sens des champs varie selon le service.")
            st.json(audit.get("usage") or {})
        with st.expander("Résultats limma et couverture documentaire"):
            rows = [{
                "Gène original": gene["original_gene_label"],
                "Symbole approuvé": gene["approved_symbol"],
                "Effet (échelle log2 normalisée)": gene["effect_log2_scale"],
                "FDR": gene["fdr_bh_treat"],
                "Documentation": STATUS_NAMES.get(gene["literature_status"], gene["literature_status"]),
            } for gene in context["genes"]]
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
            st.caption("Classement par amplitude absolue parmi les gènes retenus, et non par causalité ou intérêt thérapeutique. Les effets ne sont pas des rapports exacts d'expression brute.")
        a, b = st.columns(2)
        a.download_button("Télécharger la réponse", result["markdown"], file_name=f"reponse_{context['comparison_id']}_{index}.md", mime="text/markdown", key=f"md_{index}")
        b.download_button("Télécharger l’audit", json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False), file_name=f"audit_{context['comparison_id']}_{index}.json", mime="application/json", key=f"audit_{index}")


def main():
    st.set_page_config(page_title="Cancer Evidence Explorer", page_icon="🧬", layout="wide")
    st.title("Cancer Evidence Explorer")
    st.write("Explore les différences d’expression des gènes et les publications associées à un cancer.")
    st.caption("Prototype de recherche · Résultats TCGA et résumés scientifiques · Aucun conseil médical")

    try:
        if "builder" not in st.session_state:
            with st.spinner("Chargement des résultats et de la documentation…"):
                st.session_state.builder = RAGContextBuilder(BASE_DIR / "rag_data")
        builder = st.session_state.builder
    except (FileNotFoundError, ValueError, KeyError) as exc:
        st.error(f"Chargement impossible : {exc}")
        st.info("Place AppRAG.py à côté des trois scripts RAG et du dossier rag_data contenant les résultats annotés et literature_all_cancers.")
        st.stop()

    if "entries" not in st.session_state:
        st.session_state.entries = []
    with st.sidebar:
        st.header("Ton exploration")
        cancer = st.selectbox("Cancer étudié", builder.cancers, format_func=lambda x: CANCER_NAMES.get(x, x), key="cancer")
        mode = st.radio("Gènes à explorer", ["Les plus différentiellement exprimés", "Des gènes précis"])
        top_genes, direction, gene_symbols = 5, "both", None
        if mode == "Les plus différentiellement exprimés":
            direction = DIRECTIONS[st.selectbox("Direction d’expression", list(DIRECTIONS))]
            top_genes = st.slider("Nombre de gènes", 1, 10, 5)
        else:
            labels = sorted(builder.results.loc[builder.results.comparison_id.eq(cancer), "gene"].unique())
            gene_symbols = st.multiselect("Gènes retenus par limma (libellés originaux)", labels)
            st.caption("Sélectionne entre 1 et 10 gènes. Les anciens noms sont conservés pour la traçabilité.")
        topic = st.selectbox("Thème documentaire", list(TOPICS), format_func=lambda x: TOPICS[x][0], key="topic")
        st.divider()
        st.subheader("Modèle de langage")
        provider = st.selectbox("Fournisseur", list(PROVIDERS), key="llm_provider")
        config = PROVIDERS[provider]
        base_url = None
        token_parameter = None
        key_scope = provider
        if provider == "Compatible OpenAI":
            base_url = st.text_input("Adresse de base de ton API", placeholder="https://ton-service.example/v1", key="custom_base_url").strip()
            st.caption("Service de confiance : ta clé et le contexte seront envoyés à cette adresse. Compatible Chat Completions uniquement. Pour un serveur local, HTTP localhost est accepté.")
            key_scope += sha256(base_url.encode()).hexdigest()[:16]
            token_parameter = st.selectbox("Paramètre de limite de sortie", ["max_tokens", "max_completion_tokens"])
        model = st.text_input(
            "Identifiant exact du modèle",
            value=MODEL if provider == "Gemini" else "",
            placeholder=config["hint"], key=f"model_{provider}",
        ).strip()
        st.caption("Copie l'identifiant depuis ton fournisseur. La disponibilité et les tarifs dépendent de ton compte.")
        # A provider/endpoint change never reuses another service's key.
        api_key = st.text_input("Clé API du fournisseur choisi", type="password", key=f"api_key_{key_scope}").strip()
        if not api_key and provider == "Compatible OpenAI":
            st.caption("Laisse vide uniquement si ton serveur fonctionne sans authentification.")
        st.caption("La clé reste en mémoire de session ; elle n'est pas enregistrée dans les audits.")
        with st.expander("Volume du contexte et de la réponse"):
            articles_per_gene = st.slider("Résumés maximum par gène", 1, 3, 2)
            max_output_tokens = st.select_slider("Plafond de génération (tokens)", [1024, 2048, 4096, 8192, 16384], value=8192)
            st.caption("Ce plafond ne limite pas les tokens d'entrée et peut inclure du raisonnement selon le modèle. Réduire le nombre de gènes ou de résumés réduit le contexte.")
        if st.button("Effacer les réponses affichées"):
            st.session_state.entries = []
        st.caption("Les fichiers enregistrés sur le disque sont conservés.")

    st.subheader(CANCER_NAMES.get(cancer, cancer))
    st.write("Choisis une suggestion ou écris ta question sur la sélection à gauche.")
    c1, c2, c3 = st.columns(3)
    c1.button("Comprendre ces gènes", on_click=suggest_question, args=("explain",), width="stretch")
    c2.button("Explorer les pistes", on_click=suggest_question, args=("targets",), width="stretch")
    c3.button("Comprendre les limites", on_click=suggest_question, args=("limits",), width="stretch")
    with st.form("question_form"):
        question = st.text_area("Ta question", key="question", height=100, max_chars=4000, placeholder="Que sait-on du rôle de ces gènes dans ce cancer ?")
        submitted = st.form_submit_button("Explorer avec les sources", type="primary")
    st.caption("Chaque question est indépendante. Les réponses précédentes restent affichées mais ne sont pas envoyées au modèle. La recherche utilise des mots-clés et le thème choisi ; elle ne couvre pas toute la littérature.")

    if submitted:
        if not question.strip():
            st.warning("Écris une question avant de lancer la recherche.")
        elif not model:
            st.warning("Indique l’identifiant du modèle dans la barre latérale.")
        elif not api_key and provider != "Compatible OpenAI":
            st.warning("Ajoute la clé du fournisseur sélectionné dans la barre latérale.")
        elif mode == "Des gènes précis" and not 1 <= len(gene_symbols) <= 10:
            st.warning("Sélectionne entre 1 et 10 gènes.")
        else:
            try:
                if provider == "Compatible OpenAI":
                    validated_base_url(base_url or "")
                with st.spinner("Sélection des résultats et recherche des résumés…"):
                    context = builder.build(cancer=cancer, direction=direction, top_genes=top_genes, articles_per_gene=articles_per_gene, gene_symbols=gene_symbols, literature_query=TOPICS[topic][1] + " " + question.strip())
                    context["user_question"] = question.strip()
                    context["retrieval_topic"] = TOPICS[topic][0]
                if not context["genes"]:
                    st.info("Aucun gène retenu dans cette sélection. Aucun appel API effectué.")
                else:
                    if not context["sources"]:
                        st.info("Aucun résumé disponible : la réponse devra se limiter aux résultats et à cette absence de documentation.")
                    folder = BASE_DIR / "rag_data" / "contexts"
                    folder.mkdir(parents=True, exist_ok=True)
                    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
                    (folder / f"context_{cancer}_{stamp}.json").write_text(json.dumps(context, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
                    with st.spinner("Rédaction de la réponse sourcée…"):
                        result = answer_question(context, question.strip(), api_key=api_key, provider=provider, model=model, base_url=base_url,
                                                 max_output_tokens=max_output_tokens, token_parameter=token_parameter,
                                                 output_dir=BASE_DIR / "rag_data" / "answers")
                    if result is None:
                        st.error("Aucune réponse renvoyée par le fournisseur.")
                    else:
                        st.session_state.entries.append({"question": question.strip(), "context": context, "result": result})
            except LLMError as exc:
                st.error(str(exc))
            except ValueError as exc:
                st.error(f"Sélection ou données à vérifier : {exc}")
            except Exception:
                # Do not render SDK exception details: they can contain request data.
                st.error("La génération ou l’enregistrement a échoué. Vérifie la connexion et les droits d’écriture, puis réessaie. Aucun nouvel essai automatique.")

    entries = [(i, entry) for i, entry in enumerate(st.session_state.entries) if entry["context"]["comparison_id"] == cancer]
    for index, entry in reversed(entries):
        st.divider()
        render_answer(entry, index)
    if not entries:
        st.info("Les réponses pour ce cancer apparaîtront ici. Les tableaux détaillés et les résumés resteront repliés pour faciliter la lecture.")


if __name__ == "__main__":
    main()
