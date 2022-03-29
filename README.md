# Recrusive Partitioning Algorithm
Modifications made to the work of Rafael Durbano Lobato for use in Vention's Cartesian and Cobot palletizers.

## Requirements
- **emcc**: Emscripten compiler, which can be downloaded form a package manager of your choice (see: https://emscripten.org/docs/tools_reference/emcc.html)

## Usage
- Run `make`
- `dist/output.js` will be your output file
- In JavaScript, you can use the file like so:
```js
import { Module } from '<path_to_dist>/dist/output';

...

const packFunc = Module.cwrap('pack', 'string', ['number', 'number', 'number', 'number']);
const jsonStr = packFunc(palletLength, palletWidth, boxLength, boxWidth);

try {
  // pointList contains a list of points where your boxes can be placed.
  // Each point will be the top right corner of the box, so you may need
  // to do some math on the consumption end to whichever point you require.
  const pointList = JSON.parse(jsonStr);
  return pointList;
}
catch (e) {
  console.error(e);
  return [];
}
```

## Modifications Made
We require that the Recrusive Partitioning Algorithm be used in a JavaScript frontend application. To do this, we build it as a WebAssembly module using the emscripten compiler. Via emscripten, we expose a method called `pack`, which takes the `palletLength`, `palletWidth`, `boxLength`, and `boxWidth` as paramters and returns a json string containing the optimized box locations. We additionally output these positiosn to the `solution.mp`, which was the original implementation.