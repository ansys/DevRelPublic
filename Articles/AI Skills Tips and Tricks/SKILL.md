---
name: blog-post-editorial-review
description: Review and prepare Markdown blog posts before publication. Use for proofreading, sentence-case headings, code-fence language validation, heading hierarchy checks, and first-use acronym expansion while preserving technical meaning and Markdown.
metadata:
  version: "1.2.0"
---

# Blog post editorial review

## Inputs and outputs

Require one Markdown article. Never overwrite the source. In apply mode, create `{name}-reviewed.md` and `{name}-review-report.md`. In review-only mode, create only the report.

## Workflow

1. Preserve YAML front matter, Markdown structure, links, images, HTML, tables, blockquotes, lists, inline code, fenced code, URLs, filenames, identifiers, and quotations.
2. Correct only high-confidence spelling, grammar, punctuation, duplicated words, and typographical errors. Preserve technical meaning and authorial voice. Report ambiguity instead of guessing.
3. Convert headings to sentence case while preserving proper nouns, approved product names, acronyms, technology names, identifiers, commands, and filenames.
4. Validate fenced-code language identifiers. Preserve correct identifiers. Correct or add an identifier only when syntax or context establishes it confidently. Use conventional values such as `python`, `javascript`, `typescript`, `json`, `yaml`, `bash`, `powershell`, `html`, `css`, `sql`, `cpp`, `csharp`, `java`, `xml`, and `text`. Use `text` for logs, console output, and plain text. Leave uncertain blocks unlabelled and report them. Never editorially alter executable code.
5. Validate heading hierarchy. Permit no more than one H1. Do not allow downward skips such as H2 to H4. Correct levels only when intent is unambiguous. Otherwise preserve and report. If multiple H1 headings exist, retain the first as the title and demote later H1 headings only when the intended parent is clear.
6. Expand acronyms on first meaningful prose use as `Full term (ACRONYM)`, then use the acronym. Do not expand universally familiar terms such as HTTP, HTML, URL, API, and PDF; product names; or occurrences in code, inline code, filenames, paths, URLs, or quotations. Never guess. Report uncertain expansions under **Terminology requiring confirmation**. If first use is in a heading, prefer expansion in the first subsequent prose occurrence.
7. Verify valid Markdown, preserved source content, no more than one H1, coherent heading levels, sentence-case headings, deliberate code-fence identifiers, consistent acronym handling, and a complete report.

## Review report

Include exact counts for spelling and grammar corrections, heading-case changes, code-language changes, hierarchy corrections, acronym expansions, and manual-review items. For each change, provide its location, original text, revision, and reason. Explicitly state the H1 count and unresolved code blocks or terminology.

## Change policy

Apply only high-confidence changes that preserve meaning. Report but do not apply ambiguous wording, uncertain technical claims, unknown acronym expansions, uncertain code languages, unclear hierarchy, substantive voice changes, or changes inside quotations and executable code.

## Regression testing

Use `evals/evals.json` and its fixture files. Run each case in a clean agent session. Compare deterministic results with `expected/`. Run `python scripts/validate_output.py <reviewed-files>`. Treat altered code, URLs, YAML values, filenames, quotations, technical meaning, silently guessed uncertainties, multiple H1 headings, or skipped heading levels as failures.
