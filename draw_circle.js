function drawCircle(diameter) {
  const canvas = document.getElementById('circleCanvas');
  const ctx = canvas.getContext('2d');

  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = 'white';
  ctx.fillRect(0, 0, canvas.width, canvas.height);

  const radius = diameter / 2;
  const centerX = canvas.width / 2;
  const centerY = canvas.height / 2;

  ctx.beginPath();
  ctx.arc(centerX, centerY, radius, 0, Math.PI * 2);
  ctx.strokeStyle = 'blue';
  ctx.lineWidth = 2;
  ctx.stroke();
}

function handleDraw() {
  const input = prompt('הזן קוטר לעיגול:');

  if (input === null) {
    return;
  }

  const diameter = Number(input);

  if (!Number.isFinite(diameter) || diameter <= 0) {
    alert('נא הזן מספר חיובי.');
    return;
  }

  drawCircle(diameter);
}

window.addEventListener('load', () => {
  const button = document.getElementById('drawBtn');
  button.addEventListener('click', handleDraw);
});
