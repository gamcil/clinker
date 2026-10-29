import { similarityCsv, serializeProject, downloadText } from './project.mjs?v=403a581c8212cbeb';

const themePreferenceKey = 'clinker-interface-theme';
export function readThemePreference() {
  try {
    const saved = localStorage.getItem(themePreferenceKey);
    return ["auto", "light", "dark"].includes(saved) ? saved : "auto";
  } catch { return "auto"; }
}

// Common host controls. Analysis and file opening are supplied by each entry.
export function mountViewerControls({ viewer, editor, getTrees, getSimilarity, onEditorOpen = () => {}, onError }) {
  const toolbarMenus = [...document.querySelectorAll('.menu')];
  const helpButton = document.querySelector('#plot-help');
  const helpDialog = document.querySelector('#help-dialog');
  const editPlotButton = document.querySelector('#edit-plot');
  const appearanceLabel = document.querySelector('#appearance-label');
  const themeModeButtons = [...document.querySelectorAll('[data-theme-mode]')];
  const systemDarkMode = matchMedia('(prefers-color-scheme: dark)');
  let themePreference = readThemePreference();

  function resolvedTheme(mode) {
    return mode === "auto" ? (systemDarkMode.matches ? "dark" : "light") : mode;
  }
  function applyInterfaceTheme(mode) {
    const resolved = resolvedTheme(mode);
    document.documentElement.dataset.theme = resolved;
    const label = mode === "auto" ? "System" : mode[0].toLocaleUpperCase() + mode.slice(1);
    appearanceLabel.title = `Interface theme: ${label}`;
    appearanceLabel.setAttribute("aria-label", `Interface theme: ${label}`);
    for (const button of themeModeButtons) {
      button.setAttribute("aria-pressed", String(button.dataset.themeMode === mode));
    }
  }

  for (const menu of toolbarMenus) {
    menu.addEventListener("toggle", () => {
      if (!menu.open) return;
      toolbarMenus.forEach(other => {
        if (other !== menu) other.open = false;
      });
    });
    menu.addEventListener("click", event => {
      if (!event.target.closest("button:not(:disabled)")) return;
      menu.open = false;
      queueMicrotask(() => {
        if (menu.contains(document.activeElement)) menu.querySelector("summary").focus();
      });
    });
    menu.addEventListener("keydown", event => {
      if (event.key === "Tab") {
        setTimeout(() => {
          if (!menu.contains(document.activeElement)) menu.open = false;
        }, 0);
        return;
      }
      if (!["ArrowDown", "ArrowUp", "Home", "End"].includes(event.key)) return;
      const summary = menu.querySelector("summary");
      if (!menu.open && event.target !== summary) return;
      if (!menu.open && ["Home", "End"].includes(event.key)) return;
      event.preventDefault();
      menu.open = true;
      const items = [...menu.querySelectorAll("button:not(:disabled), input:not(:disabled)")]
        .filter(item => item.getClientRects().length);
      if (!items.length) return;
      const index = items.indexOf(document.activeElement);
      const next = event.key === "Home" ? 0 : event.key === "End" ? items.length - 1
        : event.key === "ArrowDown" ? (index + 1) % items.length
          : index < 0 ? items.length - 1 : (index - 1 + items.length) % items.length;
      items[next].focus();
    });
  }

  document.addEventListener("pointerdown", event => {
    toolbarMenus.forEach(menu => {
      if (!menu.contains(event.target)) menu.open = false;
    });
  });

  helpButton.addEventListener("click", () => {
    toolbarMenus.forEach(menu => { menu.open = false; });
    helpDialog.showModal();
  });
  // Keep plot shortcuts from changing the figure while the modal has focus.
  helpDialog.addEventListener("keydown", event => event.stopPropagation());
  helpDialog.addEventListener("close", () => helpButton.focus());
  helpDialog.addEventListener("click", event => {
    const topic = event.target.closest("[data-help-section]");
    if (topic) {
      helpDialog.querySelector(`#${topic.dataset.helpSection}`).scrollIntoView({ block: "start" });
      return;
    }
    if (event.target !== helpDialog) return;
    const bounds = helpDialog.getBoundingClientRect();
    if (event.clientX < bounds.left || event.clientX > bounds.right
        || event.clientY < bounds.top || event.clientY > bounds.bottom) helpDialog.close();
  });


    document.addEventListener('keydown', event => {
      if (event.key !== 'Escape') return;
      const openMenu = toolbarMenus.find(menu => menu.open);
      if (!openMenu) return;
      event.preventDefault();
      openMenu.open = false;
      openMenu.querySelector('summary').focus();
    });
    const syncEditorToggle = () => {
      const open = viewer.editorOpen;
      if (open) onEditorOpen();
      editPlotButton.setAttribute('aria-expanded', String(open));
      editPlotButton.textContent = open ? 'Close editor' : 'Edit figure';
    };
    editor.addEventListener('clinker-editor-open', syncEditorToggle);
    editor.addEventListener('clinker-editor-close', syncEditorToggle);
    editPlotButton.addEventListener('click', () => {
      if (viewer.editorOpen) viewer.closeEditor();
      else { onEditorOpen(); viewer.openEditor(); }
    });
    syncEditorToggle();
  document.querySelector("#download-matrix").addEventListener("click", () => {
    const latestSimilarity = getSimilarity();
    if (!latestSimilarity) return;
    try {
      downloadText("clinker-cluster-similarity.csv", similarityCsv(latestSimilarity, { data: viewer.chart.data(), state: viewer.chart.state() }), "text/csv;charset=utf-8");
    } catch (error) {
      onError(`CSV export failed: ${error.message || String(error)}`, { tone: "error" });
    }
  });

  document.querySelector("#download-svg").addEventListener("click", () => {
    try {
      downloadText("clinker.svg", viewer.exportSvg(), "image/svg+xml;charset=utf-8");
    } catch (error) {
      onError(`SVG export failed: ${error.message || String(error)}`, { tone: "error" });
    }
  });

  for (const button of themeModeButtons) {
    button.addEventListener("click", () => {
      themePreference = button.dataset.themeMode;
      viewer.setAppearance(themePreference);
      applyInterfaceTheme(themePreference);
      try { localStorage.setItem(themePreferenceKey, themePreference); } catch {}
    });
  }
  systemDarkMode.addEventListener?.("change", () => {
    if (themePreference === "auto") applyInterfaceTheme("auto");
  });

  document.querySelector("#download-project").addEventListener("click", () => {
    if (!viewer.analysis) return;
    try {
      const project = viewer.project(getTrees());
      downloadText("clinker-project.json", serializeProject(project), "application/json;charset=utf-8");
    } catch (error) {
      onError(`Project export failed: ${error.message || String(error)}`, { tone: "error" });
    }
  });
  applyInterfaceTheme(themePreference);
}
