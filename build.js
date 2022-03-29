const fs = require('fs');

fs.readFile('./dist/output.js', 'utf8', function(err, data) {
  data = data.replace("var Module = typeof Module !== 'undefined' ? Module : {};", "export var Module = typeof Module !== 'undefined' ? Module : {};");

  fs.writeFile('./dist/output.js', data, 'utf8', function(err) {
    if (err) return console.log(err);
  });
});