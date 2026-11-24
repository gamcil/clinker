import { b as uiAppearanceStyles, h as applyUiAppearance } from './search-field-8b97b3e8.mjs';

// Nested controls inherit public overrides and the enclosing editor's explicit
// appearance. Each panel also works independently with complete system defaults.
const editorThemeStyles = uiAppearanceStyles(":host") + `
:host { background:var(--cm-surface,var(--_cm-surface)); }
`;

function applyEditorTheme(element, chart) {
  applyUiAppearance(element, chart ? chart.uiAppearance() : "auto");
}

export { applyEditorTheme as a, editorThemeStyles as e };
