export function rangeShortcutKeys(platform = '') {
  if (/^(Mac|iPhone|iPad|iPod)/i.test(platform)) return { add: ['⌘'], remove: ['⌥'] };
  if (/^(Win|Linux)/i.test(platform)) return { add: ['Ctrl'], remove: ['Alt'] };
  return { add: ['Ctrl', '⌘'], remove: ['Alt', '⌥'] };
}

export function formatShortcutKeys(keys) {
  return keys.map(key => `<kbd>${key}</kbd>`).join('/');
}
