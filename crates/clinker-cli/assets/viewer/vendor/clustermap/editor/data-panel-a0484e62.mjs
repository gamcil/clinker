import { b as defineDataTable, c as defineEditableField, i as installDataEditorStyle, a as installStyle, d as defineActionMenu, e as defineSelectionBar, f as defineAnnotationEditor, R as RowSelection, D as DATA_EDITOR_CSS, g as DATA_ROW_HEIGHT } from './annotation-editor-3c63735d.mjs';
import { r as rgb$1, i as color, e as colourHex, c as createTemplate, m as mountComponent, d as updateSelectOptions, f as defineSearchField, g as createActionRunner, a as controlStyles, u as upgradeProperty } from './search-field-595834d9.mjs';
import { e as editorThemeStyles, a as applyEditorTheme } from './theme-8bd47f2f.mjs';

/** Shared activation and horizontal keyboard navigation for option buttons. */
function bindOptionButtons(group, choose) {
  group.addEventListener("click", (event) => choose(event.target.closest("button")));
  group.addEventListener("keydown", (event) => {
    if (!["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
    const buttons = [...group.querySelectorAll("button:not(:disabled)")];
    if (!buttons.length) return;
    const index = buttons.indexOf(event.target);
    const next = event.key === "Home" ? 0 : event.key === "End" ? buttons.length - 1
      : (index + (event.key === "ArrowLeft" ? -1 : 1) + buttons.length) % buttons.length;
    event.preventDefault(); buttons[next].focus(); choose(buttons[next]);
  });
}

function ascending(a, b) {
  return a == null || b == null ? NaN : a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}

function descending(a, b) {
  return a == null || b == null ? NaN
    : b < a ? -1
    : b > a ? 1
    : b >= a ? 0
    : NaN;
}

function bisector(f) {
  let compare1, compare2, delta;

  // If an accessor is specified, promote it to a comparator. In this case we
  // can test whether the search value is (self-) comparable. We can’t do this
  // for a comparator (except for specific, known comparators) because we can’t
  // tell if the comparator is symmetric, and an asymmetric comparator can’t be
  // used to test whether a single value is comparable.
  if (f.length !== 2) {
    compare1 = ascending;
    compare2 = (d, x) => ascending(f(d), x);
    delta = (d, x) => f(d) - x;
  } else {
    compare1 = f === ascending || f === descending ? f : zero$1;
    compare2 = f;
    delta = f;
  }

  function left(a, x, lo = 0, hi = a.length) {
    if (lo < hi) {
      if (compare1(x, x) !== 0) return hi;
      do {
        const mid = (lo + hi) >>> 1;
        if (compare2(a[mid], x) < 0) lo = mid + 1;
        else hi = mid;
      } while (lo < hi);
    }
    return lo;
  }

  function right(a, x, lo = 0, hi = a.length) {
    if (lo < hi) {
      if (compare1(x, x) !== 0) return hi;
      do {
        const mid = (lo + hi) >>> 1;
        if (compare2(a[mid], x) <= 0) lo = mid + 1;
        else hi = mid;
      } while (lo < hi);
    }
    return lo;
  }

  function center(a, x, lo = 0, hi = a.length) {
    const i = left(a, x, lo, hi - 1);
    return i > lo && delta(a[i - 1], x) > -delta(a[i], x) ? i - 1 : i;
  }

  return {left, center, right};
}

function zero$1() {
  return 0;
}

function number$1(x) {
  return x === null ? NaN : +x;
}

const ascendingBisect = bisector(ascending);
const bisectRight = ascendingBisect.right;
bisector(number$1).center;
var bisect = bisectRight;

class InternMap extends Map {
  constructor(entries, key = keyof) {
    super();
    Object.defineProperties(this, {_intern: {value: new Map()}, _key: {value: key}});
    if (entries != null) for (const [key, value] of entries) this.set(key, value);
  }
  get(key) {
    return super.get(intern_get(this, key));
  }
  has(key) {
    return super.has(intern_get(this, key));
  }
  set(key, value) {
    return super.set(intern_set(this, key), value);
  }
  delete(key) {
    return super.delete(intern_delete(this, key));
  }
}

function intern_get({_intern, _key}, value) {
  const key = _key(value);
  return _intern.has(key) ? _intern.get(key) : value;
}

function intern_set({_intern, _key}, value) {
  const key = _key(value);
  if (_intern.has(key)) return _intern.get(key);
  _intern.set(key, value);
  return value;
}

function intern_delete({_intern, _key}, value) {
  const key = _key(value);
  if (_intern.has(key)) {
    value = _intern.get(key);
    _intern.delete(key);
  }
  return value;
}

function keyof(value) {
  return value !== null && typeof value === "object" ? value.valueOf() : value;
}

const e10 = Math.sqrt(50),
    e5 = Math.sqrt(10),
    e2 = Math.sqrt(2);

function tickSpec(start, stop, count) {
  const step = (stop - start) / Math.max(0, count),
      power = Math.floor(Math.log10(step)),
      error = step / Math.pow(10, power),
      factor = error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1;
  let i1, i2, inc;
  if (power < 0) {
    inc = Math.pow(10, -power) / factor;
    i1 = Math.round(start * inc);
    i2 = Math.round(stop * inc);
    if (i1 / inc < start) ++i1;
    if (i2 / inc > stop) --i2;
    inc = -inc;
  } else {
    inc = Math.pow(10, power) * factor;
    i1 = Math.round(start / inc);
    i2 = Math.round(stop / inc);
    if (i1 * inc < start) ++i1;
    if (i2 * inc > stop) --i2;
  }
  if (i2 < i1 && 0.5 <= count && count < 2) return tickSpec(start, stop, count * 2);
  return [i1, i2, inc];
}

function ticks(start, stop, count) {
  stop = +stop, start = +start, count = +count;
  if (!(count > 0)) return [];
  if (start === stop) return [start];
  const reverse = stop < start, [i1, i2, inc] = reverse ? tickSpec(stop, start, count) : tickSpec(start, stop, count);
  if (!(i2 >= i1)) return [];
  const n = i2 - i1 + 1, ticks = new Array(n);
  if (reverse) {
    if (inc < 0) for (let i = 0; i < n; ++i) ticks[i] = (i2 - i) / -inc;
    else for (let i = 0; i < n; ++i) ticks[i] = (i2 - i) * inc;
  } else {
    if (inc < 0) for (let i = 0; i < n; ++i) ticks[i] = (i1 + i) / -inc;
    else for (let i = 0; i < n; ++i) ticks[i] = (i1 + i) * inc;
  }
  return ticks;
}

function tickIncrement(start, stop, count) {
  stop = +stop, start = +start, count = +count;
  return tickSpec(start, stop, count)[2];
}

function tickStep(start, stop, count) {
  stop = +stop, start = +start, count = +count;
  const reverse = stop < start, inc = reverse ? tickIncrement(stop, start, count) : tickIncrement(start, stop, count);
  return (reverse ? -1 : 1) * (inc < 0 ? 1 / -inc : inc);
}

function initRange(domain, range) {
  switch (arguments.length) {
    case 0: break;
    case 1: this.range(domain); break;
    default: this.range(range).domain(domain); break;
  }
  return this;
}

const implicit = Symbol("implicit");

function ordinal() {
  var index = new InternMap(),
      domain = [],
      range = [],
      unknown = implicit;

  function scale(d) {
    let i = index.get(d);
    if (i === undefined) {
      if (unknown !== implicit) return unknown;
      index.set(d, i = domain.push(d) - 1);
    }
    return range[i % range.length];
  }

  scale.domain = function(_) {
    if (!arguments.length) return domain.slice();
    domain = [], index = new InternMap();
    for (const value of _) {
      if (index.has(value)) continue;
      index.set(value, domain.push(value) - 1);
    }
    return scale;
  };

  scale.range = function(_) {
    return arguments.length ? (range = Array.from(_), scale) : range.slice();
  };

  scale.unknown = function(_) {
    return arguments.length ? (unknown = _, scale) : unknown;
  };

  scale.copy = function() {
    return ordinal(domain, range).unknown(unknown);
  };

  initRange.apply(scale, arguments);

  return scale;
}

var constant = x => () => x;

function linear$1(a, d) {
  return function(t) {
    return a + t * d;
  };
}

function exponential(a, b, y) {
  return a = Math.pow(a, y), b = Math.pow(b, y) - a, y = 1 / y, function(t) {
    return Math.pow(a + t * b, y);
  };
}

function gamma(y) {
  return (y = +y) === 1 ? nogamma : function(a, b) {
    return b - a ? exponential(a, b, y) : constant(isNaN(a) ? b : a);
  };
}

function nogamma(a, b) {
  var d = b - a;
  return d ? linear$1(a, d) : constant(isNaN(a) ? b : a);
}

var rgb = (function rgbGamma(y) {
  var color = gamma(y);

  function rgb(start, end) {
    var r = color((start = rgb$1(start)).r, (end = rgb$1(end)).r),
        g = color(start.g, end.g),
        b = color(start.b, end.b),
        opacity = nogamma(start.opacity, end.opacity);
    return function(t) {
      start.r = r(t);
      start.g = g(t);
      start.b = b(t);
      start.opacity = opacity(t);
      return start + "";
    };
  }

  rgb.gamma = rgbGamma;

  return rgb;
})(1);

function numberArray(a, b) {
  if (!b) b = [];
  var n = a ? Math.min(b.length, a.length) : 0,
      c = b.slice(),
      i;
  return function(t) {
    for (i = 0; i < n; ++i) c[i] = a[i] * (1 - t) + b[i] * t;
    return c;
  };
}

function isNumberArray(x) {
  return ArrayBuffer.isView(x) && !(x instanceof DataView);
}

function genericArray(a, b) {
  var nb = b ? b.length : 0,
      na = a ? Math.min(nb, a.length) : 0,
      x = new Array(na),
      c = new Array(nb),
      i;

  for (i = 0; i < na; ++i) x[i] = interpolate(a[i], b[i]);
  for (; i < nb; ++i) c[i] = b[i];

  return function(t) {
    for (i = 0; i < na; ++i) c[i] = x[i](t);
    return c;
  };
}

function date(a, b) {
  var d = new Date;
  return a = +a, b = +b, function(t) {
    return d.setTime(a * (1 - t) + b * t), d;
  };
}

function interpolateNumber(a, b) {
  return a = +a, b = +b, function(t) {
    return a * (1 - t) + b * t;
  };
}

function object(a, b) {
  var i = {},
      c = {},
      k;

  if (a === null || typeof a !== "object") a = {};
  if (b === null || typeof b !== "object") b = {};

  for (k in b) {
    if (k in a) {
      i[k] = interpolate(a[k], b[k]);
    } else {
      c[k] = b[k];
    }
  }

  return function(t) {
    for (k in i) c[k] = i[k](t);
    return c;
  };
}

var reA = /[-+]?(?:\d+\.?\d*|\.?\d+)(?:[eE][-+]?\d+)?/g,
    reB = new RegExp(reA.source, "g");

function zero(b) {
  return function() {
    return b;
  };
}

function one(b) {
  return function(t) {
    return b(t) + "";
  };
}

function string(a, b) {
  var bi = reA.lastIndex = reB.lastIndex = 0, // scan index for next number in b
      am, // current match in a
      bm, // current match in b
      bs, // string preceding current number in b, if any
      i = -1, // index in s
      s = [], // string constants and placeholders
      q = []; // number interpolators

  // Coerce inputs to strings.
  a = a + "", b = b + "";

  // Interpolate pairs of numbers in a & b.
  while ((am = reA.exec(a))
      && (bm = reB.exec(b))) {
    if ((bs = bm.index) > bi) { // a string precedes the next number in b
      bs = b.slice(bi, bs);
      if (s[i]) s[i] += bs; // coalesce with previous string
      else s[++i] = bs;
    }
    if ((am = am[0]) === (bm = bm[0])) { // numbers in a & b match
      if (s[i]) s[i] += bm; // coalesce with previous string
      else s[++i] = bm;
    } else { // interpolate non-matching numbers
      s[++i] = null;
      q.push({i: i, x: interpolateNumber(am, bm)});
    }
    bi = reB.lastIndex;
  }

  // Add remains of b.
  if (bi < b.length) {
    bs = b.slice(bi);
    if (s[i]) s[i] += bs; // coalesce with previous string
    else s[++i] = bs;
  }

  // Special optimization for only a single match.
  // Otherwise, interpolate each of the numbers and rejoin the string.
  return s.length < 2 ? (q[0]
      ? one(q[0].x)
      : zero(b))
      : (b = q.length, function(t) {
          for (var i = 0, o; i < b; ++i) s[(o = q[i]).i] = o.x(t);
          return s.join("");
        });
}

function interpolate(a, b) {
  var t = typeof b, c;
  return b == null || t === "boolean" ? constant(b)
      : (t === "number" ? interpolateNumber
      : t === "string" ? ((c = color(b)) ? (b = c, rgb) : string)
      : b instanceof color ? rgb
      : b instanceof Date ? date
      : isNumberArray(b) ? numberArray
      : Array.isArray(b) ? genericArray
      : typeof b.valueOf !== "function" && typeof b.toString !== "function" || isNaN(b) ? object
      : interpolateNumber)(a, b);
}

function interpolateRound(a, b) {
  return a = +a, b = +b, function(t) {
    return Math.round(a * (1 - t) + b * t);
  };
}

function constants(x) {
  return function() {
    return x;
  };
}

function number(x) {
  return +x;
}

var unit = [0, 1];

function identity$1(x) {
  return x;
}

function normalize(a, b) {
  return (b -= (a = +a))
      ? function(x) { return (x - a) / b; }
      : constants(isNaN(b) ? NaN : 0.5);
}

function clamper(a, b) {
  var t;
  if (a > b) t = a, a = b, b = t;
  return function(x) { return Math.max(a, Math.min(b, x)); };
}

// normalize(a, b)(x) takes a domain value x in [a,b] and returns the corresponding parameter t in [0,1].
// interpolate(a, b)(t) takes a parameter t in [0,1] and returns the corresponding range value x in [a,b].
function bimap(domain, range, interpolate) {
  var d0 = domain[0], d1 = domain[1], r0 = range[0], r1 = range[1];
  if (d1 < d0) d0 = normalize(d1, d0), r0 = interpolate(r1, r0);
  else d0 = normalize(d0, d1), r0 = interpolate(r0, r1);
  return function(x) { return r0(d0(x)); };
}

function polymap(domain, range, interpolate) {
  var j = Math.min(domain.length, range.length) - 1,
      d = new Array(j),
      r = new Array(j),
      i = -1;

  // Reverse descending domains.
  if (domain[j] < domain[0]) {
    domain = domain.slice().reverse();
    range = range.slice().reverse();
  }

  while (++i < j) {
    d[i] = normalize(domain[i], domain[i + 1]);
    r[i] = interpolate(range[i], range[i + 1]);
  }

  return function(x) {
    var i = bisect(domain, x, 1, j) - 1;
    return r[i](d[i](x));
  };
}

function copy(source, target) {
  return target
      .domain(source.domain())
      .range(source.range())
      .interpolate(source.interpolate())
      .clamp(source.clamp())
      .unknown(source.unknown());
}

function transformer() {
  var domain = unit,
      range = unit,
      interpolate$1 = interpolate,
      transform,
      untransform,
      unknown,
      clamp = identity$1,
      piecewise,
      output,
      input;

  function rescale() {
    var n = Math.min(domain.length, range.length);
    if (clamp !== identity$1) clamp = clamper(domain[0], domain[n - 1]);
    piecewise = n > 2 ? polymap : bimap;
    output = input = null;
    return scale;
  }

  function scale(x) {
    return x == null || isNaN(x = +x) ? unknown : (output || (output = piecewise(domain.map(transform), range, interpolate$1)))(transform(clamp(x)));
  }

  scale.invert = function(y) {
    return clamp(untransform((input || (input = piecewise(range, domain.map(transform), interpolateNumber)))(y)));
  };

  scale.domain = function(_) {
    return arguments.length ? (domain = Array.from(_, number), rescale()) : domain.slice();
  };

  scale.range = function(_) {
    return arguments.length ? (range = Array.from(_), rescale()) : range.slice();
  };

  scale.rangeRound = function(_) {
    return range = Array.from(_), interpolate$1 = interpolateRound, rescale();
  };

  scale.clamp = function(_) {
    return arguments.length ? (clamp = _ ? true : identity$1, rescale()) : clamp !== identity$1;
  };

  scale.interpolate = function(_) {
    return arguments.length ? (interpolate$1 = _, rescale()) : interpolate$1;
  };

  scale.unknown = function(_) {
    return arguments.length ? (unknown = _, scale) : unknown;
  };

  return function(t, u) {
    transform = t, untransform = u;
    return rescale();
  };
}

function continuous() {
  return transformer()(identity$1, identity$1);
}

function formatDecimal(x) {
  return Math.abs(x = Math.round(x)) >= 1e21
      ? x.toLocaleString("en").replace(/,/g, "")
      : x.toString(10);
}

// Computes the decimal coefficient and exponent of the specified number x with
// significant digits p, where x is positive and p is in [1, 21] or undefined.
// For example, formatDecimalParts(1.23) returns ["123", 0].
function formatDecimalParts(x, p) {
  if (!isFinite(x) || x === 0) return null; // NaN, ±Infinity, ±0
  var i = (x = p ? x.toExponential(p - 1) : x.toExponential()).indexOf("e"), coefficient = x.slice(0, i);

  // The string returned by toExponential either has the form \d\.\d+e[-+]\d+
  // (e.g., 1.2e+3) or the form \de[-+]\d+ (e.g., 1e+3).
  return [
    coefficient.length > 1 ? coefficient[0] + coefficient.slice(2) : coefficient,
    +x.slice(i + 1)
  ];
}

function exponent(x) {
  return x = formatDecimalParts(Math.abs(x)), x ? x[1] : NaN;
}

function formatGroup(grouping, thousands) {
  return function(value, width) {
    var i = value.length,
        t = [],
        j = 0,
        g = grouping[0],
        length = 0;

    while (i > 0 && g > 0) {
      if (length + g + 1 > width) g = Math.max(1, width - length);
      t.push(value.substring(i -= g, i + g));
      if ((length += g + 1) > width) break;
      g = grouping[j = (j + 1) % grouping.length];
    }

    return t.reverse().join(thousands);
  };
}

function formatNumerals(numerals) {
  return function(value) {
    return value.replace(/[0-9]/g, function(i) {
      return numerals[+i];
    });
  };
}

// [[fill]align][sign][symbol][0][width][,][.precision][~][type]
var re = /^(?:(.)?([<>=^]))?([+\-( ])?([$#])?(0)?(\d+)?(,)?(\.\d+)?(~)?([a-z%])?$/i;

function formatSpecifier(specifier) {
  if (!(match = re.exec(specifier))) throw new Error("invalid format: " + specifier);
  var match;
  return new FormatSpecifier({
    fill: match[1],
    align: match[2],
    sign: match[3],
    symbol: match[4],
    zero: match[5],
    width: match[6],
    comma: match[7],
    precision: match[8] && match[8].slice(1),
    trim: match[9],
    type: match[10]
  });
}

formatSpecifier.prototype = FormatSpecifier.prototype; // instanceof

function FormatSpecifier(specifier) {
  this.fill = specifier.fill === undefined ? " " : specifier.fill + "";
  this.align = specifier.align === undefined ? ">" : specifier.align + "";
  this.sign = specifier.sign === undefined ? "-" : specifier.sign + "";
  this.symbol = specifier.symbol === undefined ? "" : specifier.symbol + "";
  this.zero = !!specifier.zero;
  this.width = specifier.width === undefined ? undefined : +specifier.width;
  this.comma = !!specifier.comma;
  this.precision = specifier.precision === undefined ? undefined : +specifier.precision;
  this.trim = !!specifier.trim;
  this.type = specifier.type === undefined ? "" : specifier.type + "";
}

FormatSpecifier.prototype.toString = function() {
  return this.fill
      + this.align
      + this.sign
      + this.symbol
      + (this.zero ? "0" : "")
      + (this.width === undefined ? "" : Math.max(1, this.width | 0))
      + (this.comma ? "," : "")
      + (this.precision === undefined ? "" : "." + Math.max(0, this.precision | 0))
      + (this.trim ? "~" : "")
      + this.type;
};

// Trims insignificant zeros, e.g., replaces 1.2000k with 1.2k.
function formatTrim(s) {
  out: for (var n = s.length, i = 1, i0 = -1, i1; i < n; ++i) {
    switch (s[i]) {
      case ".": i0 = i1 = i; break;
      case "0": if (i0 === 0) i0 = i; i1 = i; break;
      default: if (!+s[i]) break out; if (i0 > 0) i0 = 0; break;
    }
  }
  return i0 > 0 ? s.slice(0, i0) + s.slice(i1 + 1) : s;
}

var prefixExponent;

function formatPrefixAuto(x, p) {
  var d = formatDecimalParts(x, p);
  if (!d) return prefixExponent = undefined, x.toPrecision(p);
  var coefficient = d[0],
      exponent = d[1],
      i = exponent - (prefixExponent = Math.max(-8, Math.min(8, Math.floor(exponent / 3))) * 3) + 1,
      n = coefficient.length;
  return i === n ? coefficient
      : i > n ? coefficient + new Array(i - n + 1).join("0")
      : i > 0 ? coefficient.slice(0, i) + "." + coefficient.slice(i)
      : "0." + new Array(1 - i).join("0") + formatDecimalParts(x, Math.max(0, p + i - 1))[0]; // less than 1y!
}

function formatRounded(x, p) {
  var d = formatDecimalParts(x, p);
  if (!d) return x + "";
  var coefficient = d[0],
      exponent = d[1];
  return exponent < 0 ? "0." + new Array(-exponent).join("0") + coefficient
      : coefficient.length > exponent + 1 ? coefficient.slice(0, exponent + 1) + "." + coefficient.slice(exponent + 1)
      : coefficient + new Array(exponent - coefficient.length + 2).join("0");
}

var formatTypes = {
  "%": (x, p) => (x * 100).toFixed(p),
  "b": (x) => Math.round(x).toString(2),
  "c": (x) => x + "",
  "d": formatDecimal,
  "e": (x, p) => x.toExponential(p),
  "f": (x, p) => x.toFixed(p),
  "g": (x, p) => x.toPrecision(p),
  "o": (x) => Math.round(x).toString(8),
  "p": (x, p) => formatRounded(x * 100, p),
  "r": formatRounded,
  "s": formatPrefixAuto,
  "X": (x) => Math.round(x).toString(16).toUpperCase(),
  "x": (x) => Math.round(x).toString(16)
};

function identity(x) {
  return x;
}

var map = Array.prototype.map,
    prefixes = ["y","z","a","f","p","n","µ","m","","k","M","G","T","P","E","Z","Y"];

function formatLocale(locale) {
  var group = locale.grouping === undefined || locale.thousands === undefined ? identity : formatGroup(map.call(locale.grouping, Number), locale.thousands + ""),
      currencyPrefix = locale.currency === undefined ? "" : locale.currency[0] + "",
      currencySuffix = locale.currency === undefined ? "" : locale.currency[1] + "",
      decimal = locale.decimal === undefined ? "." : locale.decimal + "",
      numerals = locale.numerals === undefined ? identity : formatNumerals(map.call(locale.numerals, String)),
      percent = locale.percent === undefined ? "%" : locale.percent + "",
      minus = locale.minus === undefined ? "−" : locale.minus + "",
      nan = locale.nan === undefined ? "NaN" : locale.nan + "";

  function newFormat(specifier, options) {
    specifier = formatSpecifier(specifier);

    var fill = specifier.fill,
        align = specifier.align,
        sign = specifier.sign,
        symbol = specifier.symbol,
        zero = specifier.zero,
        width = specifier.width,
        comma = specifier.comma,
        precision = specifier.precision,
        trim = specifier.trim,
        type = specifier.type;

    // The "n" type is an alias for ",g".
    if (type === "n") comma = true, type = "g";

    // The "" type, and any invalid type, is an alias for ".12~g".
    else if (!formatTypes[type]) precision === undefined && (precision = 12), trim = true, type = "g";

    // If zero fill is specified, padding goes after sign and before digits.
    if (zero || (fill === "0" && align === "=")) zero = true, fill = "0", align = "=";

    // Compute the prefix and suffix.
    // For SI-prefix, the suffix is lazily computed.
    var prefix = (options && options.prefix !== undefined ? options.prefix : "") + (symbol === "$" ? currencyPrefix : symbol === "#" && /[boxX]/.test(type) ? "0" + type.toLowerCase() : ""),
        suffix = (symbol === "$" ? currencySuffix : /[%p]/.test(type) ? percent : "") + (options && options.suffix !== undefined ? options.suffix : "");

    // What format function should we use?
    // Is this an integer type?
    // Can this type generate exponential notation?
    var formatType = formatTypes[type],
        maybeSuffix = /[defgprs%]/.test(type);

    // Set the default precision if not specified,
    // or clamp the specified precision to the supported range.
    // For significant precision, it must be in [1, 21].
    // For fixed precision, it must be in [0, 20].
    precision = precision === undefined ? 6
        : /[gprs]/.test(type) ? Math.max(1, Math.min(21, precision))
        : Math.max(0, Math.min(20, precision));

    function format(value) {
      var valuePrefix = prefix,
          valueSuffix = suffix,
          i, n, c;

      if (type === "c") {
        valueSuffix = formatType(value) + valueSuffix;
        value = "";
      } else {
        value = +value;

        // Determine the sign. -0 is not less than 0, but 1 / -0 is!
        var valueNegative = value < 0 || 1 / value < 0;

        // Perform the initial formatting.
        value = isNaN(value) ? nan : formatType(Math.abs(value), precision);

        // Trim insignificant zeros.
        if (trim) value = formatTrim(value);

        // If a negative value rounds to zero after formatting, and no explicit positive sign is requested, hide the sign.
        if (valueNegative && +value === 0 && sign !== "+") valueNegative = false;

        // Compute the prefix and suffix.
        valuePrefix = (valueNegative ? (sign === "(" ? sign : minus) : sign === "-" || sign === "(" ? "" : sign) + valuePrefix;
        valueSuffix = (type === "s" && !isNaN(value) && prefixExponent !== undefined ? prefixes[8 + prefixExponent / 3] : "") + valueSuffix + (valueNegative && sign === "(" ? ")" : "");

        // Break the formatted value into the integer “value” part that can be
        // grouped, and fractional or exponential “suffix” part that is not.
        if (maybeSuffix) {
          i = -1, n = value.length;
          while (++i < n) {
            if (c = value.charCodeAt(i), 48 > c || c > 57) {
              valueSuffix = (c === 46 ? decimal + value.slice(i + 1) : value.slice(i)) + valueSuffix;
              value = value.slice(0, i);
              break;
            }
          }
        }
      }

      // If the fill character is not "0", grouping is applied before padding.
      if (comma && !zero) value = group(value, Infinity);

      // Compute the padding.
      var length = valuePrefix.length + value.length + valueSuffix.length,
          padding = length < width ? new Array(width - length + 1).join(fill) : "";

      // If the fill character is "0", grouping is applied after padding.
      if (comma && zero) value = group(padding + value, padding.length ? width - valueSuffix.length : Infinity), padding = "";

      // Reconstruct the final output based on the desired alignment.
      switch (align) {
        case "<": value = valuePrefix + value + valueSuffix + padding; break;
        case "=": value = valuePrefix + padding + value + valueSuffix; break;
        case "^": value = padding.slice(0, length = padding.length >> 1) + valuePrefix + value + valueSuffix + padding.slice(length); break;
        default: value = padding + valuePrefix + value + valueSuffix; break;
      }

      return numerals(value);
    }

    format.toString = function() {
      return specifier + "";
    };

    return format;
  }

  function formatPrefix(specifier, value) {
    var e = Math.max(-8, Math.min(8, Math.floor(exponent(value) / 3))) * 3,
        k = Math.pow(10, -e),
        f = newFormat((specifier = formatSpecifier(specifier), specifier.type = "f", specifier), {suffix: prefixes[8 + e / 3]});
    return function(value) {
      return f(k * value);
    };
  }

  return {
    format: newFormat,
    formatPrefix: formatPrefix
  };
}

var locale;
var format;
var formatPrefix;

defaultLocale({
  thousands: ",",
  grouping: [3],
  currency: ["$", ""]
});

function defaultLocale(definition) {
  locale = formatLocale(definition);
  format = locale.format;
  formatPrefix = locale.formatPrefix;
  return locale;
}

function precisionFixed(step) {
  return Math.max(0, -exponent(Math.abs(step)));
}

function precisionPrefix(step, value) {
  return Math.max(0, Math.max(-8, Math.min(8, Math.floor(exponent(value) / 3))) * 3 - exponent(Math.abs(step)));
}

function precisionRound(step, max) {
  step = Math.abs(step), max = Math.abs(max) - step;
  return Math.max(0, exponent(max) - exponent(step)) + 1;
}

function tickFormat(start, stop, count, specifier) {
  var step = tickStep(start, stop, count),
      precision;
  specifier = formatSpecifier(specifier == null ? ",f" : specifier);
  switch (specifier.type) {
    case "s": {
      var value = Math.max(Math.abs(start), Math.abs(stop));
      if (specifier.precision == null && !isNaN(precision = precisionPrefix(step, value))) specifier.precision = precision;
      return formatPrefix(specifier, value);
    }
    case "":
    case "e":
    case "g":
    case "p":
    case "r": {
      if (specifier.precision == null && !isNaN(precision = precisionRound(step, Math.max(Math.abs(start), Math.abs(stop))))) specifier.precision = precision - (specifier.type === "e");
      break;
    }
    case "f":
    case "%": {
      if (specifier.precision == null && !isNaN(precision = precisionFixed(step))) specifier.precision = precision - (specifier.type === "%") * 2;
      break;
    }
  }
  return format(specifier);
}

function linearish(scale) {
  var domain = scale.domain;

  scale.ticks = function(count) {
    var d = domain();
    return ticks(d[0], d[d.length - 1], count == null ? 10 : count);
  };

  scale.tickFormat = function(count, specifier) {
    var d = domain();
    return tickFormat(d[0], d[d.length - 1], count == null ? 10 : count, specifier);
  };

  scale.nice = function(count) {
    if (count == null) count = 10;

    var d = domain();
    var i0 = 0;
    var i1 = d.length - 1;
    var start = d[i0];
    var stop = d[i1];
    var prestep;
    var step;
    var maxIter = 10;

    if (stop < start) {
      step = start, start = stop, stop = step;
      step = i0, i0 = i1, i1 = step;
    }

    while (maxIter-- > 0) {
      step = tickIncrement(start, stop, count);
      if (step === prestep) {
        d[i0] = start;
        d[i1] = stop;
        return domain(d);
      } else if (step > 0) {
        start = Math.floor(start / step) * step;
        stop = Math.ceil(stop / step) * step;
      } else if (step < 0) {
        start = Math.ceil(start * step) / step;
        stop = Math.floor(stop * step) / step;
      } else {
        break;
      }
      prestep = step;
    }

    return scale;
  };

  return scale;
}

function linear() {
  var scale = continuous();

  scale.copy = function() {
    return copy(scale, linear());
  };

  initRange.apply(scale, arguments);

  return linearish(scale);
}

Object.freeze({
  show: true, fontSize: 10, fontWeight: 600, borderWidth: 3, borderOpacity: 1,
});

const sequential = ["#f7fbff", "#08306b"];
const categorical = ["#4477aa", "#ee6677", "#228833", "#ccbb44", "#66ccee", "#aa3377", "#bbbbbb"];
Object.freeze({ cellWidth: 32, cellHeight: 20, gap: 24,
  titleFontSize: 10, columnFontSize: 10, labelFontSize: 10, headingRotation: -60, missingColour: "#eeeeee",
  cellGap: 0, cellBorderWidth: 0, cellBorderColour: "#ffffff", width: 120, colour: "#267c9e" });
const ANNOTATION_PALETTES = {
  Blue: sequential,
  "Blue–red": ["#2166ac", "#f7f7f7", "#b2182b"],
  Viridis: ["#440154", "#21918c", "#fde725"],
  Qualitative: categorical,
};

// Values are replaced at data boundaries and shared by presentation edits.
// Retain only their inferred domain, not a second flattened copy of the matrix.
const domainCache = new WeakMap();
function inferredDomain(values, type, isBar) {
  let cached = domainCache.get(values);
  if (!cached) { cached = {}; domainCache.set(values, cached); }
  const key = type === "categorical" ? "categories" : isBar ? "bars" : "numeric";
  if (cached[key]) return cached[key];
  const categories = new Set();
  let min = Infinity, max = -Infinity, count = 0;
  for (const row of Object.values(values)) for (const value of row) {
    if (value === null) continue;
    count++;
    if (type === "categorical") categories.add(value);
    else { min = Math.min(min, value); max = Math.max(max, value); }
  }
  if (isBar) { min = Math.min(0, min); max = Math.max(0, max); }
  return cached[key] = type === "categorical" ? [...categories]
    : !count || isBar && min === max ? [0, 1]
    : min === max ? [min - 0.5, max + 0.5] : [min, max];
}

/** Resolve scales identically for layout and the appearance editor. */
function resolveAnnotationScale(source) {
  const scale = source.scale ?? {}, type = scale.type ?? "continuous";
  const domain = scale.domain ?? inferredDomain(source.values, type, source.type === "bar");
  const colours = scale.colours ?? (type === "continuous" ? sequential : categorical);
  const colour = type === "categorical"
    ? ordinal(domain, colours).unknown(source.missingColour ?? "#eeeeee")
    : linear().domain(colours.map((_, i) => domain[0] + (domain[1] - domain[0]) * i / (colours.length - 1))).range(colours).clamp(true);
  return { type, domain, colours, colour };
}
Object.freeze([]);

const rowTemplates = createTemplate("<template id=\"groups\">\n  <div class=\"cm-data-row cm-editor__row cm-editor__row--groups\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\" data-select></span>\n    <button type=\"button\" class=\"cm-data-icon-button cm-editor__reorder\" draggable=\"true\" data-reorder title=\"Drag to reorder; Arrow Up or Arrow Down moves one place\">⠿</button>\n    <clinker-editable-field class=\"cm-data-cell\" data-type=\"groups\" data-field=\"label\" label=\"Group name\" placeholder=\"Group name\"></clinker-editable-field>\n    <clinker-editable-field class=\"cm-data-cell\" data-type=\"groups\" data-field=\"subtitle\" label=\"Legend label\" placeholder=\"Add subtitle…\"></clinker-editable-field>\n    <span class=\"cm-data-cell\" data-count></span>\n    <button type=\"button\" class=\"cm-data-icon-button\" data-visible>\n      <svg aria-hidden=\"true\" viewBox=\"0 0 24 24\">\n        <path d=\"M2.5 12s3.5-6 9.5-6 9.5 6 9.5 6-3.5 6-9.5 6-9.5-6Z\"></path>\n        <circle cx=\"12\" cy=\"12\" r=\"2.75\"></circle><path data-slash d=\"m4 4 16 16\"></path>\n      </svg>\n    </button>\n    <span class=\"cm-data-cell\"><input type=\"color\" data-colour aria-label=\"Colour\"></span>\n  </div>\n</template>\n\n<template id=\"links\">\n  <div class=\"cm-data-row cm-editor__row cm-editor__row--links\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\" data-select></span>\n    <span class=\"cm-data-cell\" data-text=\"source\"></span>\n    <span class=\"cm-data-cell\" data-text=\"target\"></span>\n    <clinker-editable-field class=\"cm-data-cell\" data-type=\"links\" data-field=\"identity\" label=\"Identity\" type=\"number\" min=\"0\" max=\"100\" step=\"any\" required></clinker-editable-field>\n    <span class=\"cm-data-cell\" data-text=\"groups\"></span>\n    <clinker-editable-field class=\"cm-data-cell\" data-type=\"links\" data-field=\"label\" label=\"Link label\" placeholder=\"Add label…\"></clinker-editable-field>\n    <button type=\"button\" class=\"cm-data-icon-button\" data-visible>\n      <svg aria-hidden=\"true\" viewBox=\"0 0 24 24\">\n        <path d=\"M2.5 12s3.5-6 9.5-6 9.5 6 9.5 6-3.5 6-9.5 6-9.5-6Z\"></path>\n        <circle cx=\"12\" cy=\"12\" r=\"2.75\"></circle><path data-slash d=\"m4 4 16 16\"></path>\n      </svg>\n    </button>\n    <span class=\"cm-data-cell\"><input type=\"color\" data-colour aria-label=\"Colour\"></span>\n  </div>\n</template>\n\n<template id=\"genes\">\n  <div class=\"cm-data-row cm-editor__row cm-editor__row--genes\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\" data-select></span>\n    <clinker-editable-field class=\"cm-data-cell\" data-type=\"genes\" data-field=\"label\" label=\"Gene label\" placeholder=\"Gene label\"></clinker-editable-field>\n    <span class=\"cm-data-cell\" data-text=\"locus\"></span>\n    <span class=\"cm-data-cell\" data-text=\"cluster\"></span>\n    <span class=\"cm-data-cell\" data-text=\"group\"></span>\n    <span class=\"cm-data-cell\"><input type=\"color\" data-colour aria-label=\"Colour\"></span>\n  </div>\n</template>\n\n<template id=\"gene\">\n  <div class=\"cm-data-row cm-editor__tree-row cm-editor__tree-gene\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\" data-select></span>\n    <span aria-hidden=\"true\"></span>\n    <clinker-editable-field class=\"cm-data-cell\" data-type=\"genes\" data-field=\"label\" label=\"Gene label\" placeholder=\"Gene label\"></clinker-editable-field>\n    <span class=\"cm-data-cell\" data-text=\"group\"></span>\n    <span class=\"cm-data-cell\"><input type=\"color\" data-colour aria-label=\"Colour\"></span>\n  </div>\n</template>\n\n<template id=\"cluster\">\n  <div class=\"cm-data-row cm-editor__tree-row cm-editor__tree-parent\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\"></span>\n    <button type=\"button\" class=\"cm-data-icon-button cm-editor__tree-toggle\" data-tree-toggle></button>\n    <span class=\"cm-editor__tree-name-with-count\">\n      <clinker-editable-field class=\"cm-data-cell\" data-type=\"clusters\" data-field=\"name\" label=\"Cluster name\" placeholder=\"Cluster name\" compact required></clinker-editable-field>\n      <span class=\"cm-editor__tree-count\"></span>\n      <button type=\"button\" class=\"cm-data-icon-button\" data-remove-record aria-label=\"Remove cluster from plot\" title=\"Remove cluster from plot\">×</button>\n    </span>\n  </div>\n</template>\n\n<template id=\"locus\">\n  <div class=\"cm-data-row cm-editor__tree-row cm-editor__tree-parent\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\"></span>\n    <button type=\"button\" class=\"cm-data-icon-button cm-editor__tree-toggle\" data-tree-toggle></button>\n    <span class=\"cm-editor__tree-name-with-count\">\n      <clinker-editable-field class=\"cm-data-cell\" data-type=\"loci\" data-field=\"name\" label=\"Locus name\" placeholder=\"Locus name\" compact required></clinker-editable-field>\n      <span class=\"cm-editor__tree-count\"></span>\n      <button type=\"button\" class=\"cm-data-icon-button\" data-remove-record aria-label=\"Remove locus from plot\" title=\"Remove locus from plot\">×</button>\n    </span>\n  </div>\n</template>\n\n<template id=\"cluster-locus\">\n  <div class=\"cm-data-row cm-editor__tree-row cm-editor__tree-parent\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\"></span>\n    <button type=\"button\" class=\"cm-data-icon-button cm-editor__tree-toggle\" data-tree-toggle></button>\n    <div class=\"cm-editor__tree-combined\">\n      <clinker-editable-field class=\"cm-data-cell\" data-type=\"clusters\" data-field=\"name\" label=\"Cluster name\" placeholder=\"Cluster name\" compact required></clinker-editable-field>\n      <span class=\"cm-editor__tree-separator\">·</span>\n      <clinker-editable-field class=\"cm-data-cell\" data-type=\"loci\" data-field=\"name\" label=\"Locus name\" placeholder=\"Locus name\" compact required></clinker-editable-field>\n      <span class=\"cm-editor__tree-count\"></span>\n      <button type=\"button\" class=\"cm-data-icon-button\" data-remove-record aria-label=\"Remove cluster from plot\" title=\"Remove cluster from plot\">×</button>\n    </div>\n  </div>\n</template>\n\n<template id=\"gene-group\">\n  <div class=\"cm-data-row cm-editor__tree-row cm-editor__tree-parent\">\n    <span class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\"></span>\n    <button type=\"button\" class=\"cm-data-icon-button cm-editor__tree-toggle\" data-tree-toggle></button>\n    <span class=\"cm-data-cell cm-editor__tree-group-label\"><span data-text=\"label\"></span> <span class=\"cm-editor__tree-count\"></span></span>\n  </div>\n</template>\n");
const definitions = new WeakMap();

/** Stable rows assembled from HTML templates; fields own their editing behavior. */
function createEditorRow(document, kind, type) {
  if (!definitions.has(document)) definitions.set(document, rowTemplates(document));
  const template = definitions.get(document).getElementById(type || kind);
  const node = document.importNode(template.content, true).firstElementChild;
  const parent = type && type !== "gene";
  const check = node.querySelector('[type="checkbox"]'), swatch = node.querySelector('[type="color"]');
  const visible = node.querySelector("[data-visible]"), arrow = node.querySelector("[data-tree-toggle]");
  const count = node.querySelector("[data-count], .cm-editor__tree-count");
  const reorder = node.querySelector("[data-reorder]"), fields = [...node.querySelectorAll("clinker-editable-field")];
  const texts = [...node.querySelectorAll("[data-text]")];
  return { node, update(row, selection) {
    const selected = parent ? row.selectIds.length > 0 && row.selectIds.every((id) => selection.has(id)) : selection.has(row.id);
    node.classList.toggle("is-selected", selected);
    node.style.setProperty("--depth", row.depth || 0);
    check.checked = selected;
    check.setAttribute("aria-label", parent ? `Select ${row.count.toLocaleString()} genes in ${row.label}` : `Select ${row.label || row.source}`);
    if (parent) {
      check.indeterminate = !selected && row.selectIds.some((id) => selection.has(id));
      arrow.textContent = row.open ? "▾" : "▸";
      arrow.setAttribute("aria-label", row.open ? "Collapse" : "Expand"); arrow.title = row.open ? "Collapse" : "Expand";
    }
    for (const field of fields) {
      const record = field.dataset.type === "clusters" ? row.cluster : field.dataset.type === "loci" ? row.locus : null;
      field.dataset.id = record?.uid ?? row.id;
      field.value = field.dataset.field === "name"
        ? type === "cluster-locus" && field.dataset.type === "loci" ? row.locusLabel : row.label
        : row[field.dataset.field];
    }
    for (const text of texts) text.title = text.textContent = row[text.dataset.text] ?? "";
    if (count) count.textContent = parent ? `(${row.count.toLocaleString()})` : `${row.count.toLocaleString()} genes`;
    if (swatch) swatch.value = colourHex(row.colour);
    if (visible) {
      const label = row.hidden ? "Show" : "Hide";
      visible.title = label; visible.setAttribute("aria-label", label); visible.querySelector("[data-slash]").style.display = row.hidden ? "" : "none";
    }
    if (reorder) {
      reorder.setAttribute("aria-label", `Reorder ${row.label}; use Arrow Up or Arrow Down to move it`);
    }
    return node;
  }, destroy() { for (const field of fields) field.cancelEdit(); } };
}

const [markup$3, css$3] = await Promise.all([
  "<clinker-data-table class=\"cm-data-table\" role=\"region\" aria-label=\"Editable chart data\">\n  <div slot=\"header\" data-head class=\"cm-editor__head\"></div>\n  <div data-rows></div>\n</clinker-data-table>\n",
  "clinker-record-table { display:flex; flex:1; min-width:0; min-height:0; color:inherit; font:inherit; }\nclinker-record-table[hidden] { display:none; }\nclinker-record-table .cm-editor__head, clinker-record-table .cm-editor__row { display: grid; align-items: center; gap:0; min-width:640px; padding:0; box-sizing: border-box; }\nclinker-record-table .cm-editor__head { position:sticky; top:0; z-index:2; }\nclinker-record-table .cm-editor__sort { display: flex; align-items: center; gap: 2px; min-height:24px; border:0; background:transparent; color:inherit; font-size:inherit; font-weight:inherit; text-align:left; }\nclinker-record-table .cm-editor__sort:hover { border: 0; background: transparent; color: var(--cm-accent,var(--_cm-accent,#145da8)); }\nclinker-record-table .cm-editor__sort span[aria-hidden] { color: var(--cm-muted,var(--_cm-muted,#7b8a9e)); font-size: 12px; }\nclinker-record-table .cm-editor__row { cursor: pointer; }\nclinker-record-table .cm-editor__head--groups, clinker-record-table .cm-editor__row--groups { min-width:660px; grid-template-columns: 32px 28px minmax(180px, 1fr) minmax(150px, 300px) 86px 28px 56px; }\nclinker-record-table .cm-editor__head--genes, clinker-record-table .cm-editor__row--genes { grid-template-columns: 32px minmax(160px, 1.3fr) minmax(100px, .85fr) minmax(110px, 1fr) minmax(100px, .85fr) 56px; }\nclinker-record-table .cm-editor__head--links, clinker-record-table .cm-editor__row--links { min-width:802px; grid-template-columns: 32px minmax(120px, 1fr) minmax(120px, 1fr) 76px minmax(150px, .95fr) minmax(120px, .8fr) 28px 56px; }\nclinker-record-table .cm-editor__head--tree, clinker-record-table .cm-editor__tree-row { grid-template-columns:32px 28px minmax(180px, 1fr) minmax(120px, .75fr) 56px; }\nclinker-record-table .cm-editor__tree-row { display:grid; align-items:center; min-width:640px; gap:0; padding:0; }\nclinker-record-table .cm-editor__tree-parent > :last-child { grid-column:3 / -1; margin-left:calc(var(--depth, 0) * 18px); }\nclinker-record-table .cm-editor__tree-gene > clinker-editable-field { padding-left:calc(8px + var(--depth, 0) * 18px); }\nclinker-record-table .cm-editor__tree-parent { --cm-row-background:var(--cm-subtle,var(--_cm-subtle,#f3f6fa)); font-weight: 550; }\nclinker-record-table .cm-editor__tree-toggle { --cm-icon-size:24px; font-size:14px; }\nclinker-record-table .cm-editor__tree-separator { color: var(--cm-muted,var(--_cm-muted,#96a4b5)); font-size: 12px; }\nclinker-record-table .cm-editor__tree-combined { display: flex; min-width: 0; align-items: center; gap: 4px; padding-inline:8px; }\nclinker-record-table .cm-editor__tree-combined > clinker-editable-field { padding-inline:0; }\nclinker-record-table .cm-editor__tree-combined > [data-type=\"clusters\"] { color:var(--cm-foreground,var(--_cm-foreground,#324155)); font-weight:650; }\nclinker-record-table .cm-editor__tree-combined > [data-type=\"loci\"] { color:var(--cm-muted,var(--_cm-muted,#526178)); font-weight:550; }\nclinker-record-table .cm-editor__tree-count { color: var(--cm-muted,var(--_cm-muted,#5f6d7f)); font-size: 12px; font-weight: 600; white-space: nowrap; }\nclinker-record-table .cm-editor__tree-name-with-count { display: flex; align-items: center; gap: 5px; min-width: 0; }\nclinker-record-table .cm-editor__reorder { --cm-icon-size:24px; cursor:grab; font-size:15px; letter-spacing:-3px; }\nclinker-record-table .cm-editor__drop-target { box-shadow: inset 0 2px var(--cm-accent,var(--_cm-accent,#1677ff)); }\n",
]);
const template$3 = createTemplate(markup$3);
const HTMLElementBase$3 = globalThis.HTMLElement || class {};

class RecordTable extends HTMLElementBase$3 {
  #table;
  #head;
  #list;
  #state;
  #kind = null;
  #records = new WeakMap();
  #dragged = null;
  #hoverRow = null;
  #hoveredRow = null;
  #hoverTimer = null;
  #hoverFrame = null;
  error = "";

  constructor() {
    super();
    this.addEventListener("click", (event) => this.#click(event));
    this.addEventListener("change", (event) => {
      const input = event.target, node = input.closest(".cm-data-row"), row = this.#records.get(node)?.data;
      if (!row) return;
      if (input.matches('[type="checkbox"]')) {
        const checked = row.selectIds ? row.selectIds.every((id) => this.#state.selection.has(id)) : this.#state.selection.has(row.id);
        if (checked !== input.checked) this.#select(row, input.checked, false);
      } else if (input.matches("[data-colour]")) this.#emit({ action: "update", type: this.#kind, id: row.id, changes: { colour: input.value } });
    });
    this.addEventListener("value-commit", (event) => {
      if (event.defaultPrevented) return;
      const { type, field } = event.target.dataset, row = this.#records.get(event.target.closest(".cm-data-row"))?.data;
      if (!row || !type || !field) return;
      const id = type === "clusters" ? row.cluster.uid : type === "loci" ? row.locus.uid : row.id;
      event.stopPropagation();
      const value = field === "identity" ? Number(event.detail.value) / 100 : event.detail.value.trim();
      if (!this.#emit({ action: "update", type, id, changes: { [field]: value } }, event.target)) event.preventDefault();
    });
    this.addEventListener("keydown", (event) => {
      const grip = event.target.closest("[data-reorder]");
      if (!grip || event.defaultPrevented || !["ArrowUp", "ArrowDown"].includes(event.key)) return;
      event.preventDefault();
      this.#emit({ action: "reorder", id: this.#records.get(grip.closest(".cm-data-row")).data.id, direction: event.key === "ArrowUp" ? -1 : 1 });
    });
    this.addEventListener("dragstart", (event) => {
      const grip = event.target.closest("[data-reorder]");
      if (!grip) return;
      this.#dragged = this.#records.get(grip.closest(".cm-data-row")).data.id;
      event.dataTransfer.effectAllowed = "move"; event.dataTransfer.setData("text/plain", String(this.#dragged));
    });
    this.addEventListener("dragover", (event) => {
      const target = event.target.closest(".cm-editor__row--groups");
      if (this.#dragged === null || !target || this.#records.get(target).data.id === this.#dragged) return;
      event.preventDefault(); this.#clearDropTarget(); target.classList.add("cm-editor__drop-target");
    });
    this.addEventListener("drop", (event) => {
      const target = event.target.closest(".cm-editor__row--groups");
      if (this.#dragged !== null && target) {
        event.preventDefault(); const bounds = target.getBoundingClientRect();
        this.#emit({ action: "place", id: this.#dragged, targetId: this.#records.get(target).data.id, after: event.clientY > bounds.top + bounds.height / 2 });
      }
      this.#dragged = null; this.#clearDropTarget();
    });
    this.addEventListener("dragend", () => { this.#dragged = null; this.#clearDropTarget(); });
    for (const type of ["pointerover", "pointerout"]) this.addEventListener(type, (event) => {
      const row = event.target.closest(".cm-data-row");
      if (row && !row.contains(event.relatedTarget)) this.#queueHover(type === "pointerover" ? row : null);
    });
  }

  connectedCallback() {
    const root = this.getRootNode();
    installDataEditorStyle(root);
    this.classList.add("cm-ui", "cm-data-ui");
    installStyle(root, "clinker-record-table-style", css$3);
    if (!this.#table) {
      this.append(template$3(this.ownerDocument));
      this.#table = this.querySelector("clinker-data-table"); this.#head = this.querySelector("[data-head]"); this.#list = this.querySelector("[data-rows]");
    }
    this.#render();
  }
  disconnectedCallback() { this.cancelHover(); this.#dragged = null; this.#clearDropTarget(); }
  set state(value) { this.#state = value; this.#render(); }
  reset() { this.cancelHover(); if (this.#table) this.#table.renderer = null; this.#kind = null; this.#dragged = null; this.#clearDropTarget(); }
  get scrollTop() { return this.#table?.scrollTop || 0; }
  set scrollTop(value) { if (this.#table) this.#table.scrollTop = value; }
  get scrollHeight() { return this.#table?.scrollHeight || 0; }
  scrollTo(...args) { this.#table?.scrollTo(...args); }

  #emit(detail, control) {
    this.error = "";
    const accepted = this.dispatchEvent(new CustomEvent("record-action", { detail, bubbles: true, composed: true, cancelable: true }));
    if (!accepted && control) control.setCustomValidity(this.error);
    return accepted;
  }
  #select(row, checked, range) { this.#emit({ action: "select", id: row.key ?? row.id, checked, range }); }
  #click(event) {
    if (event.target.closest("clinker-editable-field")) return;
    const node = event.target.closest(".cm-data-row"), row = this.#records.get(node)?.data;
    if (row && event.target.matches('[type="checkbox"]')) { this.#select(row, event.target.checked, event.shiftKey); return; }
    const sort = event.target.closest("[data-sort]");
    if (sort) { this.#emit({ action: "sort", key: sort.dataset.sort }); return; }
    if (!row) return;
    if (event.target.closest("[data-tree-toggle]")) this.#emit({ action: "toggle", key: row.key });
    else if (event.target.closest("[data-visible]")) this.#emit({ action: "update", type: this.#kind, id: row.id, changes: { hidden: !row.hidden } });
    else if (event.target.closest("[data-remove-record]")) {
      const cluster = row.type === "cluster" || row.type === "cluster-locus";
      this.#emit({ action: "remove", type: cluster ? "clusters" : "loci", id: row[cluster ? "cluster" : "locus"].uid });
    } else if (!event.target.closest("input, select, button")) this.#select(row, !node.querySelector('[type="checkbox"]').checked, event.shiftKey);
  }
  #render() {
    if (!this.#table || !this.#state) return;
    const { kind, geneView, rows, sort, attached = true } = this.#state;
    if (!attached) {
      this.reset(); this.#table.header = { columns: [], sort: null };
      this.#list.textContent = "Set the chart property to edit a ClusterMap chart."; return;
    }
    const tree = kind === "genes" && geneView !== "list", layout = tree ? "tree" : kind;
    const columns = tree ? [{ label: "" }, { label: "" }, { label: geneView === "groups" ? "Group / gene" : "Cluster / locus / gene" }, { label: "Group" }, { label: "Colour" }]
      : kind === "groups" ? [{ label: "" }, { label: "↕" }, { label: "Group" }, { label: "Legend label" }, { label: "Genes" }, { label: "" }, { label: "Colour" }]
      : kind === "links" ? [{ label: "" }, { label: "Source", sort: "source" }, { label: "Target", sort: "target" }, { label: "Identity", sort: "identity" }, { label: "Endpoint groups", sort: "groups" }, { label: "Link label", sort: "label" }, { label: "" }, { label: "Colour" }]
        : [{ label: "" }, { label: "Gene", sort: "label" }, { label: "Locus", sort: "locus" }, { label: "Cluster", sort: "cluster" }, { label: "Group", sort: "group" }, { label: "Colour" }];
    this.#head.className = `cm-data-head cm-editor__head cm-editor__head--${layout}`;
    this.#table.header = { columns, sort }; this.#list.className = `cm-editor__rows cm-editor__rows--${layout}`;
    if (this.#kind !== kind) {
      this.reset(); this.#kind = kind; this.#list.replaceChildren();
      const updateRow = (node, data) => { const row = this.#records.get(node); row.data = data; return row.view.update(data, this.#state.selection); };
      this.#table.renderer = { container: this.#list, key: (row) => `${kind}:${row.type || "row"}:${row.key ?? row.id}`,
        headerHeight: () => this.#head.offsetHeight,
        createRow: (data) => { const view = createEditorRow(this.ownerDocument, kind, data.type); this.#records.set(view.node, { view, data }); return updateRow(view.node, data); },
        updateRow, disposeRow: (node) => this.#records.get(node).view.destroy() };
    }
    this.#table.rows = rows;
  }

  #clearDropTarget() { for (const node of this.querySelectorAll(".cm-editor__drop-target")) node.classList.remove("cm-editor__drop-target"); }
  #queueHover(row) {
    this.#hoverRow = row;
    if (this.#hoverTimer !== null) globalThis.clearTimeout(this.#hoverTimer);
    this.#hoverTimer = globalThis.setTimeout(() => {
      this.#hoverTimer = null;
      if (this.#hoverFrame !== null) return;
      const frame = globalThis.requestAnimationFrame || ((callback) => globalThis.setTimeout(callback, 0));
      this.#hoverFrame = frame(() => {
        this.#hoverFrame = null;
        const hovered = this.#hoverRow;
        if (hovered === this.#hoveredRow) return;
        this.#hoveredRow = hovered;
        this.#emit({ action: "hover", row: hovered?.isConnected ? this.#records.get(hovered).data : null });
      });
    }, 36);
  }
  cancelHover() {
    if (this.#hoverTimer !== null) globalThis.clearTimeout(this.#hoverTimer);
    if (this.#hoverFrame !== null) (globalThis.cancelAnimationFrame || globalThis.clearTimeout)(this.#hoverFrame);
    this.#hoverTimer = this.#hoverFrame = this.#hoverRow = this.#hoveredRow = null;
  }
}

function defineRecordTable() {
  defineDataTable(); defineEditableField();
  if (!customElements.get("clinker-record-table")) customElements.define("clinker-record-table", RecordTable);
}

const emptyData = { clusters: [], groups: [], links: [] };
const labelFor = (record, fallback = "") => record?.label || record?.name || fallback;

/** Cached projections of chart data, independent of DOM and chart mutations. */
class DataPanelModel {
  #data = null;
  #index = null;
  #view = null;

  get data() { return this.#data; }
  set data(value) { this.#data = value; this.#index = this.#view = null; }
  groupId(value) { return this.#records().groupById.get(String(value))?.uid; }
  recordIds(kind) {
    return kind === "genes" ? this.#records().genes.map(({ gene }) => gene.uid)
      : (this.#data?.[kind] || []).map(({ uid }) => uid);
  }

  plotTarget(kind, ids) {
    if (kind === "links") return { links: [...ids] };
    if (kind === "groups") {
      const selected = new Set(ids);
      ids = (this.#data?.groups || []).filter((group) => selected.has(group.uid)).flatMap((group) => group.genes || []);
    }
    return { genes: [...ids] };
  }

  selectedParents(kind, selection) {
    if (!selection.size) return [];
    const clusters = this.#records().clusters;
    const parents = kind === "loci" ? clusters.flatMap(({ loci }) => loci) : clusters;
    return parents.filter(({ ids }) => ids.length && ids.every((id) => selection.has(id))).map(({ record }) => record.uid);
  }

  view(options) {
    const key = JSON.stringify([options.kind, options.geneView, options.filter, options.memberGroupId, options.sort, [...options.expanded]]);
    if (this.#view?.key === key) return this.#view;
    const branches = [], memberGroup = this.#data?.groups.find((group) => group.uid === options.memberGroupId);
    const rows = this.#computeRows({ ...options, branches, memberGroup });
    const ids = [...new Set(rows.flatMap((row) => row.selectIds || (row.type === "gene" || !row.type ? [row.id] : [])))];
    return this.#view = { key, rows, ids, branches, memberGroup, hasCollapsed: branches.some((row) => !row.open), byKey: new Map(rows.map((row) => [row.key ?? row.id, row])),
      parentCount: rows.filter((row) => row.type && row.type !== "gene").length };
  }

  #records() {
    if (this.#index) return this.#index;
    const data = this.#data || emptyData;

    const groups = new Map((data.groups || []).flatMap((group) => (group.genes || []).map((id) => [id, group])));
    const groupById = new Map((data.groups || []).map((group) => [String(group.uid), group]));
    const genes = [], clusters = [];
    const geneById = new Map();
    for (const cluster of data.clusters || []) {
      const loci = [];
      for (const locus of cluster.loci || []) {
        const ids = [];
        for (const gene of locus.genes || []) {
          geneById.set(gene.uid, gene);
          genes.push({ gene, locus, cluster });
          ids.push(gene.uid);
        }
        loci.push({ record: locus, ids });
      }
      clusters.push({ record: cluster, loci, ids: loci.flatMap(({ ids }) => ids) });
    }
    return this.#index = { groups, groupById, genes, geneById, clusters };
  }

  #computeRows({ kind, geneView, filter, memberGroupId, memberGroup, sort, expanded, branches }) {
    const data = this.#data;
    if (!data) return [];
    const needle = filter;
    const includes = (...values) => !needle || values.join(" ").toLowerCase().includes(needle);
    const cache = this.#records();
    const groups = cache.groups;
    if (kind === "groups") return data.groups.filter((group) => includes(labelFor(group, group.uid), group.subtitle)).map((group) => ({
      id: group.uid, label: labelFor(group, group.uid), subtitle: group.subtitle || "", count: group.genes?.length || 0, hidden: Boolean(group.hidden), colour: group.colour,
    }));
    if (kind === "links") {
      const genes = cache.geneById;
      return sortRows(data.links.filter((link) => includes(labelFor(genes.get(link.query.uid), link.query.uid), labelFor(genes.get(link.target.uid), link.target.uid), link.label)).map((link) => ({
        id: link.uid,
        source: labelFor(genes.get(link.query.uid), link.query.uid),
        target: labelFor(genes.get(link.target.uid), link.target.uid),
        identity: Math.round((link.identity || 0) * 1000) / 10,
        groups: `${labelFor(groups.get(link.query.uid), "Ungrouped")} / ${labelFor(groups.get(link.target.uid), "Ungrouped")}`,
        label: link.label || "",
        hidden: Boolean(link.hidden),
        colour: link.colour,
      })), sort);
    }
    memberGroupId = memberGroup?.uid ?? null;
    const memberIds = memberGroup && new Set(memberGroup.genes || []);
    const genes = cache.genes.filter(({ gene, locus, cluster }) =>
      (!memberIds || memberIds.has(gene.uid)) && includes(labelFor(gene, gene.uid), labelFor(locus, locus.uid), labelFor(cluster, cluster.uid), labelFor(groups.get(gene.uid), "Ungrouped")));
    if (geneView === "tree") return treeRows(cache.clusters, groups, genes, { needle, memberGroupId, expanded, branches });
    if (geneView === "groups") return groupTreeRows(data, groups, genes, { needle, expanded, branches });
    return sortRows(genes.map(({ gene, locus, cluster }) => geneRow(gene, locus, cluster, groups)), sort);
  }
}

function sortRows(rows, sort) {
  if (!sort) return rows;
  return rows.sort((left, right) => {
    const a = left[sort.key] ?? "";
    const b = right[sort.key] ?? "";
    const comparison = typeof a === "number" && typeof b === "number"
      ? a - b
      : String(a).localeCompare(String(b), undefined, { numeric: true, sensitivity: "base" });
    return comparison * sort.direction;
  });
}

function geneRow(gene, locus, cluster, groups, depth = null) {
  return { ...(depth === null ? {} : { type: "gene", depth }), id: gene.uid, label: labelFor(gene, gene.uid), locus: labelFor(locus, locus.uid), cluster: labelFor(cluster, cluster.uid), group: labelFor(groups.get(gene.uid), "Ungrouped"), colour: gene.colour || groups.get(gene.uid)?.colour };
}

function treeRows(clusters, groups, visible, { needle, memberGroupId, expanded, branches }) {
  const visibleIds = new Set(visible.map(({ gene }) => gene.uid));
  const rows = [];
  for (const { record: cluster, loci: members } of clusters) {
    const loci = members.map(({ record: locus }) => ({ locus, genes: locus.genes.filter((gene) => visibleIds.has(gene.uid)) }))
      .filter(({ locus, genes }) => genes.length || (!needle && memberGroupId == null && !locus.genes.length));
    if (!loci.length && (cluster.loci.length || needle || memberGroupId != null)) continue;
    if (loci.length === 1) {
      const { locus, genes } = loci[0];
      const key = `cluster-locus:${cluster.uid}:${locus.uid}`;
      const open = Boolean(needle) || expanded.has(key);
      const row = { type: "cluster-locus", key, cluster, locus, label: labelFor(cluster, cluster.uid), locusLabel: labelFor(locus, locus.uid), count: genes.length, selectIds: genes.map((gene) => gene.uid), open };
      branches.push(row); rows.push(row);
      if (open) rows.push(...genes.map((gene) => geneRow(gene, locus, cluster, groups, 1)));
      continue;
    }
    const key = `cluster:${cluster.uid}`;
    const open = Boolean(needle) || expanded.has(key);
    const row = { type: "cluster", key, cluster, label: labelFor(cluster, cluster.uid), count: loci.reduce((count, item) => count + item.genes.length, 0), selectIds: loci.flatMap(({ genes }) => genes.map((gene) => gene.uid)), open };
    branches.push(row); rows.push(row);
    for (const { locus, genes } of loci) {
      const locusKey = `locus:${locus.uid}`;
      const locusOpen = Boolean(needle) || expanded.has(locusKey);
      const row = { type: "locus", key: locusKey, parentKey: key, cluster, locus, label: labelFor(locus, locus.uid), count: genes.length, selectIds: genes.map((gene) => gene.uid), open: locusOpen, depth: 1 };
      branches.push(row);
      if (!open) continue;
      rows.push(row);
      if (locusOpen) rows.push(...genes.map((gene) => geneRow(gene, locus, cluster, groups, 2)));
    }
  }
  return rows;
}

function groupTreeRows(data, groups, visible, { needle, expanded, branches }) {
  const entries = new Map(visible.map((entry) => [entry.gene.uid, entry]));
  const buckets = [
    ...(data.groups || []).map((group) => ({ group, key: `gene-group:${group.uid}`, label: labelFor(group, group.uid), genes: (group.genes || []).map((id) => entries.get(id)).filter(Boolean) })),
    { key: "gene-group:ungrouped", label: "Ungrouped", genes: visible.filter(({ gene }) => !groups.has(gene.uid)) },
  ];
  const rows = [];
  for (const bucket of buckets) {
    if (!bucket.genes.length) continue;
    const open = Boolean(needle) || expanded.has(bucket.key);
    const row = { type: "gene-group", group: bucket.group, key: bucket.key, label: bucket.label, count: bucket.genes.length, selectIds: bucket.genes.map(({ gene }) => gene.uid), open };
    branches.push(row); rows.push(row);
    if (open) rows.push(...bucket.genes.map(({ gene, locus, cluster }) => geneRow(gene, locus, cluster, groups, 1)));
  }
  return rows;
}

const [markup$2, css$2] = await Promise.all([
  "<div role=\"group\">\n  <template><button type=\"button\"></button></template>\n</div>\n",
  ":host { display: block; color: inherit; font: inherit; }\n:host([hidden]) { display: none !important; }\ndiv {\n  display: flex;\n  padding: 3px;\n  border: 1px solid var(--cm-border,var(--_cm-border, #e0e5eb));\n  border-radius: 7px;\n  background: var(--cm-subtle,var(--_cm-subtle, #f3f6fa));\n}\nbutton {\n  min-height: var(--cm-control-height,var(--_cm-control-height,28px));\n  border: 0;\n  border-radius: var(--cm-control-radius,var(--_cm-control-radius, 6px));\n  background: transparent;\n  white-space: nowrap;\n}\nbutton[aria-pressed=true] { background: var(--cm-surface,var(--_cm-surface, #fff)); box-shadow: 0 1px 2px #1d2b3d29; font-weight: 600; }\nbutton:hover:not(:disabled) { background: var(--cm-hover,var(--_cm-hover, #e3ecf8)); }\n",
]);
const template$2 = createTemplate(markup$2);
const HTMLElementBase$2 = globalThis.HTMLElement || class {};

class ViewSwitcher extends HTMLElementBase$2 {
  static observedAttributes = ["value", "label"];
  #group;
  #template;
  #buttons = [];
  #observer;

  constructor() {
    super();
    const root = mountComponent(this, template$2, css$2);
    this.#group = root.querySelector("div"); this.#template = root.querySelector("template");
    bindOptionButtons(this.#group, (button) => this.#choose(button));
    this.#observer = new this.ownerDocument.defaultView.MutationObserver((records) => {
      if (records.some((record) => record.target !== this || record.type === "childList")) this.#renderOptions();
    });
  }

  get value() { return this.getAttribute("value") || ""; }
  set value(value) { if (this.value !== String(value)) this.setAttribute("value", String(value)); }
  attributeChangedCallback() { this.#sync(); }
  connectedCallback() {
    this.#renderOptions();
    this.#observer.observe(this, { childList: true, subtree: true, characterData: true, attributes: true, attributeFilter: ["value", "label", "disabled"] });
  }
  disconnectedCallback() { this.#observer.disconnect(); }

  #renderOptions() {
    this.#buttons = [...this.querySelectorAll("option")].map((option) => {
      const button = this.ownerDocument.importNode(this.#template.content, true).firstElementChild;
      button.value = option.value; button.textContent = option.label; button.disabled = option.disabled;
      return button;
    });
    this.#group.replaceChildren(this.#template, ...this.#buttons);
    if (!this.hasAttribute("value")) this.value = this.#buttons.find((button) => !button.disabled)?.value || "";
    this.#sync();
  }

  #sync() {
    this.#group.setAttribute("aria-label", this.getAttribute("label") || "View");
    const active = this.#buttons.find((button) => !button.disabled && button.value === this.value)
      || this.#buttons.find((button) => !button.disabled);
    for (const button of this.#buttons) {
      button.setAttribute("aria-pressed", String(button.value === this.value));
      button.tabIndex = button === active ? 0 : -1;
    }
  }

  #choose(button) {
    if (!button || button.disabled || button.value === this.value) return;
    this.value = button.value;
    this.dispatchEvent(new CustomEvent("value-change", { detail: { value: this.value }, bubbles: true, composed: true }));
  }
}

function defineViewSwitcher() {
  if (!customElements.get("clinker-view-switcher")) customElements.define("clinker-view-switcher", ViewSwitcher);
}

const [markup$1, css$1] = await Promise.all([
  "<button type=\"button\" data-command=\"focus\">Focus in plot</button>\n<clinker-action-menu label=\"Assign…\" data-menu=\"assign\">\n  <label>Move to <select data-assign-group aria-label=\"Move to group\"></select></label>\n  <button type=\"button\" data-command=\"assign\">Move</button>\n  <button type=\"button\" data-command=\"unassign\">Remove from groups</button>\n  <button type=\"button\" data-command=\"new-group-for-selection\">New group…</button>\n</clinker-action-menu>\n<clinker-action-menu label=\"Edit…\" data-menu=\"edit\">\n  <input data-batch-label placeholder=\"Set label\" aria-label=\"Set label\">\n  <button type=\"button\" data-command=\"batch-label\">Set label</button>\n  <label>Set colour <input type=\"color\" value=\"#888888\" aria-label=\"Set colour\" data-batch-colour></label>\n  <select data-batch-visibility aria-label=\"Set visibility\"><option value=\"\">Visibility…</option><option value=\"show\">Show</option><option value=\"hide\">Hide</option></select>\n  <input type=\"number\" data-batch-identity min=\"0\" max=\"100\" step=\".1\" placeholder=\"Set identity (%)\" aria-label=\"Set identity (%)\">\n  <button type=\"button\" data-command=\"batch-identity\">Set identity</button>\n</clinker-action-menu>\n<clinker-action-menu label=\"Merge…\" data-menu=\"merge\">\n  <select data-merge-target aria-label=\"Group to keep\"></select>\n  <button type=\"button\" data-command=\"merge\">Merge groups</button>\n</clinker-action-menu>\n<clinker-action-menu label=\"Group…\" data-menu=\"group\">\n  <button type=\"button\" data-command=\"view-members\">View member genes</button>\n  <label>Position <input type=\"number\" data-group-position required min=\"1\"></label>\n</clinker-action-menu>\n<clinker-action-menu label=\"Remove…\" data-menu=\"remove\">\n  <button type=\"button\" data-command=\"delete\" class=\"danger\"></button>\n  <button type=\"button\" data-command=\"remove-loci\"></button>\n  <button type=\"button\" data-command=\"remove-clusters\"></button>\n  <button type=\"button\" data-command=\"delete-groups\" title=\"Genes remain in the plot\"></button>\n</clinker-action-menu>\n",
  ":host { display: flex; align-items: center; gap: 6px; color: inherit; font: inherit; }\ninput[type=number] { width: 112px; }\n.danger { border-color: var(--cm-danger,var(--_cm-danger, #e4aaa5)); color: var(--cm-danger,var(--_cm-danger, #b3261e)); background: var(--cm-hover,var(--_cm-hover, #fffafa)); }\n",
]);
const template$1 = createTemplate(markup$1);
const groupOptions = (groups) => groups.map(({ uid, label, name }) => [uid, label || name || uid]);
const HTMLElementBase$1 = globalThis.HTMLElement || class {};

class RecordActions extends HTMLElementBase$1 {
  #context;

  constructor() {
    super();
    const root = mountComponent(this, template$1, css$1);
    root.addEventListener("click", (event) => {
      const command = event.target.closest("[data-command]")?.dataset.command;
      if (!command) return;
      const field = { assign: "assign-group", merge: "merge-target", "batch-label": "batch-label", "batch-identity": "batch-identity" }[command];
      const input = field && this.#field(field);
      if (input && !input.reportValidity()) return;
      if (input && !input.value.trim()) return;
      const changes = command === "batch-label" ? { label: input.value.trim() }
        : command === "batch-identity" ? { identity: input.valueAsNumber / 100 } : null;
      this.#emit(changes ? "update" : command, changes || input?.value);
    });
    root.addEventListener("change", (event) => {
      const input = event.target;
      if (input.matches("[data-batch-colour]")) this.#emit("update", { colour: input.value });
      else if (input.matches("[data-batch-visibility]") && input.value) { this.#emit("update", { hidden: input.value === "hide" }); input.value = ""; }
      else if (input.matches("[data-group-position]") && input.reportValidity()) this.#emit("group-position", Number(input.value) - 1);
    });
    root.addEventListener("keydown", (event) => {
      if (event.key === "Enter" && event.target.matches("[data-group-position]")) event.target.blur();
    });
  }

  connectedCallback() {
    if (this.#context === undefined) this.state = { kind: "groups", count: 0, groups: [], allGroups: [], loci: [], clusters: [] };
  }

  set state({ kind, count, groups, allGroups, loci, clusters }) {
    const context = JSON.stringify([kind, count > 0, count > 1, groups.map(({ uid }) => uid)]);
    if (this.#context !== context) {
      this.#context = context;
      for (const menu of this.shadowRoot.querySelectorAll("clinker-action-menu")) menu.open = false;
      this.#field("batch-label").value = this.#field("batch-identity").value = "";
    }
    this.hidden = !count;
    const visible = { assign: kind === "genes", edit: kind !== "genes" && count > 1, merge: groups.length > 1, group: groups.length === 1, remove: true };
    for (const menu of this.shadowRoot.querySelectorAll("[data-menu]")) menu.hidden = !visible[menu.dataset.menu];
    this.#field("batch-identity").hidden = this.#command("batch-identity").hidden = kind !== "links";
    updateSelectOptions(this.#field("assign-group"), [["", "Choose a group…"], ...groupOptions(allGroups)]);
    updateSelectOptions(this.#field("merge-target"), [["", "Choose group to keep…"], ...groupOptions(groups)]);
    const position = this.#field("group-position"); position.max = allGroups.length || 1;
    if (this.shadowRoot.activeElement !== position) position.value = allGroups.findIndex(({ uid }) => uid === groups[0]?.uid) + 1;
    for (const [command, label, show] of [["delete", `Delete selected ${kind}`, true], ["remove-loci", `Entire loci (${loci.length})`, loci.length],
      ["remove-clusters", `Entire clusters (${clusters.length})`, clusters.length], ["delete-groups", `Delete groups (${groups.length})`, kind === "genes" && groups.length]]) {
      const button = this.#command(command); button.textContent = label; button.hidden = !show;
    }
  }

  #field(name) { return this.shadowRoot.querySelector(`[data-${name}]`); }
  #command(name) { return this.shadowRoot.querySelector(`[data-command="${name}"]`); }
  #emit(command, value) { this.dispatchEvent(new CustomEvent("selection-action", { detail: { command, value }, bubbles: true, composed: true })); }
}

function defineRecordActions() {
  defineActionMenu();
  if (!customElements.get("clinker-record-actions")) customElements.define("clinker-record-actions", RecordActions);
}

const [markup, css] = await Promise.all([
  "<div class=\"cm-ui cm-data-ui cm-data-layout\" part=\"content\">\n<clinker-view-switcher data-kind-picker label=\"Data type\" value=\"groups\">\n  <option value=\"groups\">Groups</option><option value=\"genes\">Genes</option>\n  <option value=\"links\">Links</option><option value=\"annotations\">Annotations</option>\n</clinker-view-switcher>\n<div part=\"toolbar\" data-record-toolbar class=\"cm-editor__toolbar cm-data-toolbar\">\n  <clinker-view-switcher data-gene-view label=\"Gene view\" value=\"tree\" hidden>\n    <option value=\"tree\">By locus</option><option value=\"groups\">By group</option><option value=\"list\">List</option>\n  </clinker-view-switcher>\n  <button type=\"button\" data-command=\"toggle-tree\" hidden>Expand visible</button>\n  <button type=\"button\" data-command=\"new-group\">Create group…</button>\n  <clinker-search-field data-filter label=\"Filter records\" placeholder=\"Search names, labels, loci, or clusters\"></clinker-search-field>\n</div>\n<clinker-annotation-editor embedded hidden></clinker-annotation-editor>\n<div data-record-editor class=\"cm-data-panel\">\n  <div data-record-tools class=\"cm-data-tools\" hidden>\n    <div class=\"cm-editor__member-filter\" data-member-filter hidden>\n      <span>Group: <strong data-member-label></strong></span>\n      <button type=\"button\" data-command=\"clear-member-filter\" aria-label=\"Clear group member filter\" title=\"Clear group member filter\">×</button>\n    </div>\n  </div>\n  <clinker-selection-bar part=\"selection\" class=\"cm-data-meta\"><div slot=\"actions\" class=\"cm-data-actions\"><clinker-record-actions></clinker-record-actions></div></clinker-selection-bar>\n  <clinker-record-table part=\"table\" class=\"cm-editor__table\"></clinker-record-table>\n</div>\n<div data-editor-error class=\"cm-data-error\" role=\"alert\" hidden></div>\n</div>\n",
  ":host { display:flex; flex:1; min-width:0; min-height:0; }\n.cm-data-layout { display:grid; flex:1; min-width:0; min-height:0; grid-template-columns:auto minmax(0,1fr); grid-template-rows:auto minmax(0,1fr) auto; grid-template-areas:\"views tools\" \"content content\" \"error error\"; gap:8px; }\n.cm-data-layout > [data-kind-picker] { grid-area:views; align-self:center; }\n.cm-data-layout > [data-record-toolbar] { grid-area:tools; min-width:0; }\n.cm-data-layout > [data-record-editor] { grid-area:content; }\n.cm-data-layout > [data-editor-error] { grid-area:error; }\n.cm-editor__member-filter { display: flex; align-items: center; gap: 6px; color: var(--cm-muted,var(--_cm-muted,#526178)); font-size: 12px; }\n.cm-editor__member-filter button { display: grid; place-items: center; width: 24px; min-height: 24px; padding: 0; border: 0; background: transparent; color: var(--cm-muted,var(--_cm-muted,#64748b)); font-size: 17px; line-height: 1; }\n@media (max-width:700px) {\n  .cm-data-layout { grid-template-columns:minmax(0,1fr); grid-template-rows:auto auto minmax(0,1fr) auto; grid-template-areas:\"views\" \"tools\" \"content\" \"error\"; }\n}\n",
]);
const template = createTemplate(markup);

const HTMLElementBase = globalThis.HTMLElement || class {};

class ClinkerDataPanelElement extends HTMLElementBase {
  #chart = null;
  #runEdit = createActionRunner((message) => {
    const status = this.shadowRoot.querySelector("[data-editor-error]");
    if (status) { status.textContent = message; status.hidden = !message; }
    if (this.#table) this.#table.error = message;
  });
  #annotationEditor = null;
  #unsubscribe = null;
  #editing = false;
  #pendingEdit = null;
  #kind = "groups";
  #geneView = "tree";
  #selected = new RowSelection();
  #table = null;
  #actions = null;
  #selectionBar = null;
  #model = new DataPanelModel();
  #filter = "";
  #memberGroupId = null;
  #expanded = new Set();
  #sort = {
    genes: { key: "label", direction: 1 },
    links: { key: "source", direction: 1 },
  };

  constructor() {
    super();
    const root = this.attachShadow({ mode: "open" });
    const style = this.ownerDocument.createElement("style"); style.id = "cm-data-editor-style"; style.textContent = controlStyles + DATA_EDITOR_CSS + editorThemeStyles + css;
    root.append(style);
    this.#bind();
    upgradeProperty(this, "chart");
  }

  connectedCallback() {
    this.classList.add("cm-ui", "cm-data-ui");
    if (!this.#table) this.#renderShell();
    this.#subscribe(); this.#syncTheme(); this.#refreshData();
    const edit = this.#pendingEdit; this.#pendingEdit = null; edit?.();
  }

  disconnectedCallback() {
    this.#unsubscribe?.(); this.#unsubscribe = null;
    this.#table?.reset();
    this.#chart?.highlight?.([], { owner: this });
  }

  get chart() { return this.#chart; }
  set chart(value) {
    if (value === this.#chart) return;
    this.#unsubscribe?.(); this.#unsubscribe = null;
    this.#table?.cancelHover(); this.#chart?.highlight?.([], { owner: this });
    this.#chart = value || null;
    this.#subscribe(); this.#syncTheme(); this.#refreshData({ replace: true });
  }

  #syncTheme() {
    applyEditorTheme(this, this.#chart);
  }

  #subscribe() {
    if (this.#unsubscribe || !this.isConnected || !this.#chart?.on) return;
    this.#unsubscribe = this.#chart.on("change", (change) => {
      if (["data.apply", "data.replace", "state.replace"].includes(change.type)) {
        if (this.#editing) return;
        this.#refreshData({ replace: change.type === "data.replace" });
      } else if (change.type === "theme.change") this.#syncTheme();
    });
  }

  editAnnotations(panelUid) {
    if (!this.isConnected) { this.#pendingEdit = () => this.editAnnotations(panelUid); return this; }
    this.#kind = "annotations";
    this.#selected.clear();
    this.#renderData();
    if (panelUid) this.#annotationEditor?.select(panelUid);
    return this;
  }

  editSelection({ kind, ids }) {
    if (!this.isConnected) {
      const selection = { kind, ids: [...ids] };
      this.#pendingEdit = () => this.editSelection(selection); return this;
    }
    this.#kind = "genes"; this.#geneView = "tree";
    this.#filter = ""; this.#memberGroupId = null;
    const selected = new Set(ids), genes = [];
    for (const row of this.#view().branches) {
      if (!selected.has(row[kind === "clusters" ? "cluster" : "locus"]?.uid)) continue;
      this.#expanded.add(row.key);
      if (row.parentKey) this.#expanded.add(row.parentKey);
      if (row.locus) genes.push(...row.selectIds);
    }
    this.#selected.clear(); this.#selected.setMany(genes);
    this.#renderData();
    const first = this.#view().rows.findIndex((row) => row.selectIds?.some((id) => this.#selected.has(id)));
    this.shadowRoot.querySelector(".cm-editor__table").scrollTop = Math.max(0, first) * DATA_ROW_HEIGHT;
    return this;
  }

  #renderShell() {
    this.shadowRoot.append(template(this.ownerDocument));
    this.#annotationEditor = this.shadowRoot.querySelector("clinker-annotation-editor");
    this.#table = this.shadowRoot.querySelector("clinker-record-table");
    this.#table.addEventListener("record-action", (event) => this.#recordAction(event));
    this.#selectionBar = this.shadowRoot.querySelector("clinker-selection-bar");
    this.#selectionBar.addEventListener("clear-selection", () => this.#clearDataSelection());
    this.#selectionBar.addEventListener("select-all", (event) => {
      this.#selected.setMany(event.detail.ids); this.#renderData();
    });
    this.#actions = this.shadowRoot.querySelector("clinker-record-actions");
    this.#actions.addEventListener("selection-action", (event) => this.#command(event.detail.command, event.detail.value));
  }

  #bind() {
    this.shadowRoot.addEventListener("click", (event) => {
      const command = event.target.closest("[data-command]")?.dataset.command;
      if (command) this.#command(command);
    });
    this.shadowRoot.addEventListener("value-change", (event) => {
      if (event.target.matches("[data-kind-picker]")) {
        this.#kind = event.detail.value;
        if (this.#kind !== "genes") this.#memberGroupId = null;
        this.#selected.clear();
      } else if (event.target.matches("[data-gene-view]")) this.#geneView = event.detail.value;
      else return;
      this.shadowRoot.querySelector(".cm-editor__table")?.scrollTo({ top: 0 });
      this.#renderData();
    });
    this.shadowRoot.addEventListener("input", (event) => {
      if (event.target.matches("[data-filter]")) this.#setDataFilter(event.target.value);
    });
    this.shadowRoot.addEventListener("keydown", (event) => {
      if (event.defaultPrevented) return;
      if (event.key === "Escape") this.#clearDataSelection();
    });
  }

  #recordAction(event) {
    if (event.defaultPrevented) return;
    const { action, id, key, type, changes, row } = event.detail;
    if (action === "select") this.#selectRow(id, event.detail.checked, event.detail.range);
    else if (action === "sort") {
      const current = this.#sort[this.#kind];
      this.#sort[this.#kind] = { key, direction: current?.key === key ? -current.direction : 1 };
      this.#table.scrollTo({ top: 0 }); this.#renderData();
    } else if (action === "toggle") {
      if (this.#expanded.has(key)) this.#expanded.delete(key); else this.#expanded.add(key);
      this.#renderData();
    } else if (action === "update") {
      if (!this.#apply({ type: `${type}.update`, ids: [id], changes })) event.preventDefault();
    } else if (action === "remove") this.#removeRecords(type, [id]);
    else if (action === "reorder") this.#moveGroup(id, this.#model.data.groups.findIndex((group) => group.uid === id) + event.detail.direction);
    else if (action === "place") this.#placeGroup(id, event.detail.targetId, event.detail.after);
    else if (action === "hover") this.#highlight(row);
  }

  #command(command, value) {
    if (command === "remove-loci" || command === "remove-clusters") {
      const kind = command === "remove-loci" ? "loci" : "clusters";
      this.#removeRecords(kind, this.#model.selectedParents(kind, this.#selected));
      return;
    }
    if (command === "delete" || command === "delete-groups") {
      if (command === "delete-groups") this.#deleteSelected("groups", this.#selectedGroups().map((group) => group.uid));
      else this.#deleteSelected();
    } else if (command === "new-group" || command === "new-group-for-selection") {
      const label = window.prompt("Name for the new homology group");
      if (label?.trim()) {
        const geneIds = command === "new-group-for-selection" ? [...this.#selected] : [];
        this.#apply({ type: "groups.create", group: { uid: this.#nextGroupUid(), label: label.trim() }, geneIds }, geneIds.length > 0);
      }
    } else if (command === "assign") {
      if (value) this.#apply({ type: "groups.assignGenes", groupId: this.#model.groupId(value), geneIds: [...this.#selected] });
    } else if (command === "unassign") {
      this.#apply({ type: "groups.unassignGenes", geneIds: [...this.#selected] });
    } else if (command === "focus") {
      if (this.#selected.size) this.#chart?.focus?.(this.#model.plotTarget(this.#kind, this.#selected));
    } else if (command === "update") {
      this.#apply({ type: `${this.#kind}.update`, ids: [...this.#selected], changes: value });
    } else if (command === "merge") {
      const targetId = this.#model.groupId(value);
      if (targetId !== undefined) {
        const sourceIds = this.#selectedGroups().map((group) => group.uid);
        this.#apply({ type: "groups.merge", targetId, sourceIds }, true);
      }
    } else if (command === "group-position") {
      this.#moveGroup(this.#selectedGroups()[0]?.uid, value);
    } else if (command === "view-members") {
      const groupId = this.#selectedGroups()[0]?.uid;
      if (groupId !== undefined) this.#viewGroupMembers(groupId);
    } else if (command === "toggle-tree") {
      this.#toggleTree();
    } else if (command === "clear-member-filter") {
      this.#memberGroupId = null;
      this.#renderData();
    }
  }

  #apply(operation, clearSelection = false) {
    if (!this.#chart?.patch) return false;
    return this.#edit(() => this.#chart.patch([operation]), clearSelection);
  }

  #edit(action, clearSelection) {
    this.#editing = true;
    const accepted = this.#runEdit(() => { action(); if (clearSelection) this.#selected.clear(); });
    this.#editing = false;
    if (accepted) this.#refreshData();
    return accepted;
  }

  #refreshData({ replace = false } = {}) {
    const data = this.#chart?.visibleData?.() || this.#chart?.data?.() || null;
    this.#model.data = data;
    if (!data?.groups.some((group) => group.uid === this.#memberGroupId)) this.#memberGroupId = null;
    if (replace) {
      this.#annotationEditor?.reset();
      this.#table?.reset();
      this.#selected.clear(); this.#expanded.clear(); this.#memberGroupId = null;
      this.#initializeExpanded();
    } else this.#selected.prune(this.#model.recordIds(this.#kind), this.#view().byKey.keys());
    this.#renderData();
  }

  #initializeExpanded() {
    const geneCount = this.#model.recordIds("genes").length;
    if (!geneCount || geneCount > 500) return;
    this.#expandTree("tree"); this.#expandTree("groups");
  }

  #expandTree(view) {
    for (const row of this.#view({ kind: "genes", geneView: view }).branches) this.#expanded.add(row.key);
  }

  #nextGroupUid() {
    const used = new Set((this.#model.data?.groups || []).map((group) => String(group.uid)));
    let index = used.size + 1;
    while (used.has(`group-${index}`)) index += 1;
    return `group-${index}`;
  }

  #selectRow(id, checked, range) {
    const { byKey } = this.#view();
    this.#selected.select(id, checked, { range, orderedIds: [...byKey.keys()],
      idsFor: (key) => byKey.get(key).selectIds || [byKey.get(key).id] });
    this.#renderData();
  }

  #view(options = {}) {
    return this.#model.view({ kind: this.#kind, geneView: this.#geneView, filter: this.#filter,
      memberGroupId: this.#memberGroupId, sort: this.#sort[this.#kind], expanded: this.#expanded, ...options });
  }

  #setDataFilter(text) {
    this.#filter = text.trim().toLowerCase();
    this.shadowRoot.querySelector(".cm-editor__table")?.scrollTo({ top: 0 });
    this.#renderData();
  }

  #clearDataSelection() {
    if (this.#kind === "annotations") this.#annotationEditor?.clearSelection();
    else { this.#selected.clear(); this.#renderData(); }
  }

  #renderData() {
    if (!this.isConnected || !this.#table) return;
    this.#table.cancelHover(); this.#highlight();
    const annotations = this.#kind === "annotations";
    this.shadowRoot.querySelector("[data-kind-picker]").value = this.#kind;
    this.shadowRoot.querySelector("[data-gene-view]").value = this.#geneView;
    this.#annotationEditor.hidden = !annotations;
    this.shadowRoot.querySelector("[data-record-toolbar]").hidden = annotations;
    this.shadowRoot.querySelector("[data-filter]").value = this.#filter;
    this.shadowRoot.querySelector("[data-record-editor]").hidden = annotations;
    if (annotations) {
      this.#table?.reset();
      this.#annotationEditor.state = { data: this.#chart?.data?.() || { clusters: [] }, annotationApi: { ANNOTATION_PALETTES, resolveAnnotationScale }, clusters: this.#model.data?.clusters || [],
        patch: (operations) => this.#chart.patch(operations) };
      return;
    }
    const data = this.#model.data;

    if (!data) {
      this.#table.state = { attached: false };
      this.#selectionBar.state = { selection: this.#selected, ids: [], summary: "No chart attached.", noun: this.#kind.slice(0, -1) };
      this.#actions.state = { kind: this.#kind, count: 0, groups: [], allGroups: [], loci: [], clusters: [] };
      return;
    }
    const geneView = this.shadowRoot.querySelector("[data-gene-view]");
    const memberFilter = this.shadowRoot.querySelector("[data-member-filter]");
    const { rows, ids, parentCount, memberGroup, hasCollapsed } = this.#view();
    memberFilter.hidden = !memberGroup || this.#kind !== "genes";
    this.shadowRoot.querySelector("[data-record-tools]").hidden = memberFilter.hidden;
    if (memberGroup) this.shadowRoot.querySelector("[data-member-label]").textContent = memberGroup.label || memberGroup.name || memberGroup.uid;
    geneView.hidden = this.#kind !== "genes";
    this.shadowRoot.querySelector('[data-command="new-group"]').hidden = this.#kind !== "groups";
    const isTree = this.#kind === "genes" && this.#geneView !== "list";
    const treeButton = this.shadowRoot.querySelector('[data-command="toggle-tree"]');
    treeButton.hidden = !isTree || Boolean(this.#filter);
    if (isTree) {
      treeButton.textContent = hasCollapsed ? "Expand visible" : "Collapse all";
      treeButton.title = hasCollapsed ? "Expand all currently visible branches" : "Collapse every tree branch";
    }
    this.#table.state = { kind: this.#kind, geneView: this.#geneView, rows, selection: this.#selected, sort: this.#sort[this.#kind] };
    const summaryText = isTree
      ? `${parentCount.toLocaleString()} ${this.#geneView === "groups" ? "group" : "cluster/locus"} rows · ${ids.length.toLocaleString()} genes${this.#filter ? " matching search" : ""}.`
      : `${rows.length.toLocaleString()} visible ${this.#kind}${memberGroup && this.#kind === "genes" ? " in selected group" : ""}${this.#filter ? " matching search" : ""}.`;
    this.#selectionBar.state = { selection: this.#selected, ids, summary: summaryText,
      noun: this.#kind.slice(0, -1), selectSuffix: isTree ? " genes" : "" };
    const parentRemoval = this.#kind === "genes" && this.#geneView === "tree" && this.#chart?.remove;
    this.#actions.state = { kind: this.#kind, count: this.#selected.size, groups: this.#selectedGroups(), allGroups: data.groups,
      loci: parentRemoval ? this.#model.selectedParents("loci", this.#selected) : [], clusters: parentRemoval ? this.#model.selectedParents("clusters", this.#selected) : [] };
  }

  #selectedGroups() {
    if (this.#kind === "groups") return (this.#model.data?.groups || []).filter((group) => this.#selected.has(group.uid));
    if (this.#kind !== "genes" || this.#geneView !== "groups") return [];
    return this.#view().rows.filter((row) => row.type === "gene-group" && row.group?.genes?.length &&
      row.group.genes.every((id) => this.#selected.has(id))).map((row) => row.group);
  }

  #removeRecords(kind, ids) {
    if (!ids.length || !window.confirm(`Remove ${ids.length} ${kind} and all their genes from the plot? This can be undone.`)) return;
    this.#edit(() => this.#chart.remove({ [kind]: ids }), true);
  }

  #deleteSelected(kind = this.#kind, ids = [...this.#selected]) {
    if (!ids.length || !window.confirm(`Delete ${ids.length} selected ${kind}?${kind === "groups" ? " Their genes will remain in the plot." : ""}`)) return;
    this.#apply({ type: `${kind}.delete`, ids }, true);
  }

  #moveGroup(uid, index) {
    const ids = this.#model.data.groups.map((group) => group.uid), from = ids.indexOf(uid);
    if (from < 0 || !Number.isFinite(index)) return;
    const to = Math.max(0, Math.min(ids.length - 1, Math.trunc(index)));
    if (to === from) return;
    ids.splice(from, 1);
    ids.splice(to, 0, uid);
    this.#apply({ type: "groups.reorder", ids });
  }

  #viewGroupMembers(groupId) {
    this.#memberGroupId = groupId;
    this.#kind = "genes";
    this.#geneView = "list";
    this.#selected.clear();
    this.shadowRoot.querySelector(".cm-editor__table")?.scrollTo({ top: 0 });
    this.#renderData();
  }

  #placeGroup(uid, targetUid, after) {
    const ids = this.#model.data.groups.map((group) => group.uid);
    const from = ids.indexOf(uid), target = ids.indexOf(targetUid);
    if (from < 0 || target < 0 || uid === targetUid) return;
    this.#moveGroup(uid, target + Number(after) - Number(from < target));
  }

  #toggleTree() {
    if (this.#view().hasCollapsed) this.#expandTree(this.#geneView);
    else this.#expanded.clear();
    this.#renderData();
  }

  #highlight(row) {
    this.#chart?.highlight?.(this.#model.plotTarget(this.#kind, row ? row.selectIds || [row.id] : this.#selected), { owner: this });
  }

}

function defineDataPanel() {
  defineSelectionBar(); defineViewSwitcher(); defineSearchField(); defineRecordActions(); defineRecordTable(); defineAnnotationEditor();
  if (!customElements.get("clinker-data-panel")) customElements.define("clinker-data-panel", ClinkerDataPanelElement);
}

export { bindOptionButtons as b, defineDataPanel as d };
