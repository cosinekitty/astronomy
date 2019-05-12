# JavaScript examples for the browser
You can use the JavaScript version of the 
[Astronomy Engine](https://cosinekitty.github.io/astronomy)
to perform client-side astronomy calculations in a web browser.
Calculations are offloaded to the visitor's computer.
Just grab a copy of 
[astronomy.js](https://github.com/cosinekitty/astronomy/blob/master/source/js/astronomy.js)
and save it on your server. Inside your HTML code, pull in the script as usual:

```html
<script src="astronomy.js"></script>
```

All the functionality is wrapped inside an object called `Astronomy`.

Here are some example web pages using Astronomy Engine in a web browser.

---

### [Moon Phase Calculator](moonphase.html)
This example shows how to determine the Moon's current phase,
and how to predict when the next few quarter phases will occur.

---

# JavaScript examples for Node.js
The same JavaScript source file
[astronomy.js](https://github.com/cosinekitty/astronomy/blob/master/source/js/astronomy.js)
works as a Node.js module. Download the file into your project directory. 
Then in your own source file, do this:

```javascript
const Astronomy = require('astronomy.js');
```

---

# [API Reference](../../source/js/README.md)
Complete documentation for all the functions and types available
in the JavaScript version of the Astronomy Engine.
