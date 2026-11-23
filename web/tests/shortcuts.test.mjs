import test from 'node:test';
import assert from 'node:assert/strict';
import { rangeShortcutKeys } from '../src/viewer/shortcuts.mjs';

test('range shortcut labels follow the reported platform, with a fallback', () => {
  assert.deepEqual(rangeShortcutKeys('MacIntel'), { add: ['⌘'], remove: ['⌥'] });
  assert.deepEqual(rangeShortcutKeys('Win32'), { add: ['Ctrl'], remove: ['Alt'] });
  assert.deepEqual(rangeShortcutKeys('Linux x86_64'), { add: ['Ctrl'], remove: ['Alt'] });
  assert.deepEqual(rangeShortcutKeys(''), { add: ['Ctrl', '⌘'], remove: ['Alt', '⌥'] });
});
