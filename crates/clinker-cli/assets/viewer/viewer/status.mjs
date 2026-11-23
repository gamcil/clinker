export function createPlotStatus(element) {
  const messageElement = element.querySelector('[role="status"]');
  const dismiss = element.querySelector('button');
  let timer;
  const clear = () => {
    clearTimeout(timer);
    element.hidden = true;
    messageElement.textContent = '';
    delete element.dataset.tone;
  };
  dismiss.addEventListener('click', clear);
  return (message, { tone = '' } = {}) => {
    if (!message) return clear();
    if (element.dataset.tone === 'error' && tone !== 'error') return;
    clearTimeout(timer);
    element.dataset.tone = tone;
    messageElement.textContent = message;
    dismiss.hidden = tone !== 'error';
    element.hidden = false;
    if (tone !== 'error') timer = setTimeout(clear, 5000);
  };
}
