#include <math.h>
#include <string.h>
#include "ch.h"
#include "hal.h"
#include "chprintf.h"
#include "nanovna.h"

//#define __DRAW_Z__

#define SWAP(x,y) { int z=x; x = y; y = z; }

static void cell_draw_marker_info(int m, int n, int w, int h);
static void frequency_string(char *buf, size_t len, int32_t freq, char *prefix);
static void markmap_all_markers(void);

//#define GRID_COLOR 0x0863
//uint16_t grid_color = 0x1084;

/* indicate dirty cells */
static uint16_t markmap[2][8];
static uint16_t current_mappage = 0;

static int32_t fgrid = 50000000;
static int16_t grid_offset;
static int16_t grid_width;

int area_width = AREA_WIDTH_NORMAL;
int area_height = HEIGHT+1;

#define GRID_RECTANGULAR (1<<0)
#define GRID_SMITH       (1<<1)
#define GRID_ADMIT       (1<<2)
#define GRID_POLAR       (1<<3)

#define CELLWIDTH 32
#define CELLHEIGHT 32

//#define floatToInt(v) ((int)(v))
int floatToInt(float v){
	if (v < 0) return v-0.5;
	if (v > 0) return v+0.5;
	return 0;
}

/*
 * CELL_X0[27:31] cell position
 * CELL_Y0[22:26]
 * CELL_N[10:21] original order
 * CELL_X[5:9] position in the cell
 * CELL_Y[0:4]
 */
static uint32_t trace_index[TRACE_COUNT][POINT_COUNT];

#define INDEX(x, y, n) \
  ((((x)&0x03e0UL)<<22) | (((y)&0x03e0UL)<<17) | (((n)&0x0fffUL)<<10)  \
 | (((x)&0x1fUL)<<5) | ((y)&0x1fUL))

#define CELL_X(i) (int)((((i)>>5)&0x1f) | (((i)>>22)&0x03e0))
#define CELL_Y(i) (int)(((i)&0x1f) | (((i)>>17)&0x03e0))
#define CELL_N(i) (int)(((i)>>10)&0xfff)

#define CELL_X0(i) (int)(((i)>>22)&0x03e0)
#define CELL_Y0(i) (int)(((i)>>17)&0x03e0)

#define CELL_P(i, x, y) (((((x)&0x03e0UL)<<22) | (((y)&0x03e0UL)<<17)) == ((i)&0xffc00000UL))


void update_grid(void)
{
  int32_t gdigit = 100000000;
  int32_t fstart, fspan;
  int32_t grid;
  if (frequency1 > 0) {
    fstart = frequency0;
    fspan = frequency1 - frequency0;
  } else {
    fspan = -frequency1;
    fstart = frequency0 - fspan/2;
  }
  
  while (gdigit > 100) {
    grid = 5 * gdigit;
    if (fspan / grid >= 4)
      break;
    grid = 2 * gdigit;
    if (fspan / grid >= 4)
      break;
    grid = gdigit;
    if (fspan / grid >= 4)
      break;
    gdigit /= 10;
  }
  fgrid = grid;

  grid_offset = (WIDTH-1) * ((fstart % fgrid) / 100) / (fspan / 100);
  grid_width = (WIDTH-1) * (fgrid / 100) / (fspan / 1000);

  force_set_markmap();
  redraw_request |= REDRAW_FREQUENCY;
}

static int circle_inout(int x, int y, int r)
{
  int d = x*x + y*y - r*r;
  if (d <= -r)
    return 1;
  if (d > r)
    return -1;
  return 0;
}

static int polar_grid(int x, int y)
{
  int d;

  // offset to center
  x -= P_CENTER_X;
  y -= P_CENTER_Y;

  // outer circle
  d = circle_inout(x, y, P_RADIUS);
  if (d < 0) return 0;
  if (d == 0) return 1;

  // vertical and horizontal axis
  if (x == 0 || y == 0)
    return 1;

  d = circle_inout(x, y, P_RADIUS / 5);
  if (d == 0) return 1;
  if (d > 0) return 0;

  d = circle_inout(x, y, P_RADIUS * 2 / 5);
  if (d == 0) return 1;
  if (d > 0) return 0;

  // cross sloping lines
  if (x == y || x == -y)
    return 1;

  d = circle_inout(x, y, P_RADIUS * 3 / 5);
  if (d == 0) return 1;
  if (d > 0) return 0;

  d = circle_inout(x, y, P_RADIUS * 4 / 5);
  if (d == 0) return 1;
  return 0;
}

/*
 * Constant Resistance circle: (u - r/(r+1))^2 + v^2 = 1/(r+1)^2
 * Constant Reactance circle:  (u - 1)^2 + (v-1/x)^2 = 1/x^2
 */
static int smith_grid(int x, int y)
{
  int d;

  // offset to center
  x -= P_CENTER_X;
  y -= P_CENTER_Y;
  
  // outer circle
  d = circle_inout(x, y, P_RADIUS);
  if (d < 0)
    return 0;
  if (d == 0)
    return 1;
  
  // horizontal axis
  if (y == 0)
    return 1;

  // shift circle center to right origin
  x -= P_RADIUS;

  // Constant Reactance Circle: 2j : R/2 = 58
  if (circle_inout(x, y+P_RADIUS/2, P_RADIUS/2) == 0)
    return 1;
  if (circle_inout(x, y-P_RADIUS/2, P_RADIUS/2) == 0)
    return 1;

  // Constant Resistance Circle: 3 : R/4 = 29
  d = circle_inout(x+P_RADIUS/4, y, P_RADIUS/4);
  if (d > 0) return 0;
  if (d == 0) return 1;

  // Constant Reactance Circle: 1j : R = 116
  if (circle_inout(x, y+P_RADIUS, P_RADIUS) == 0)
    return 1;
  if (circle_inout(x, y-P_RADIUS, P_RADIUS) == 0)
    return 1;

  // Constant Resistance Circle: 1 : R/2 = 58
  d = circle_inout(x+P_RADIUS/2, y, P_RADIUS/2);
  if (d > 0) return 0;
  if (d == 0) return 1;

  // Constant Reactance Circle: 1/2j : R*2 = 232
  if (circle_inout(x, y+P_RADIUS*2, P_RADIUS*2) == 0)
    return 1;
  if (circle_inout(x, y-P_RADIUS*2, P_RADIUS*2) == 0)
    return 1;

  // Constant Resistance Circle: 1/3 : R*3/4 = 87
  if (circle_inout(x+P_RADIUS*3/4, y, P_RADIUS*3/4) == 0)
    return 1;
  return 0;
}

//static int smith_grid2(int x, int y, float scale)
//{
//  int d;
//
//  // offset to center
//  x -= P_CENTER_X;
//  y -= P_CENTER_Y;
//  
//  // outer circle
//  d = circle_inout(x, y, P_RADIUS);
//  if (d < 0)
//    return 0;
//  if (d == 0)
//    return 1;
//
//  // shift circle center to right origin
//  x -= P_RADIUS * scale;
//
//  // Constant Reactance Circle: 2j : R/2 = 58
//  if (circle_inout(x, y+58*scale, 58*scale) == 0)
//    return 1;
//  if (circle_inout(x, y-58*scale, 58*scale) == 0)
//    return 1;
//#if 0
//  // Constant Resistance Circle: 3 : R/4 = 29
//  d = circle_inout(x+29*scale, y, 29*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//  d = circle_inout(x-29*scale, y, 29*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//#endif
//
//  // Constant Reactance Circle: 1j : R = 116
//  if (circle_inout(x, y+116*scale, 116*scale) == 0)
//    return 1;
//  if (circle_inout(x, y-116*scale, 116*scale) == 0)
//    return 1;
//
//  // Constant Resistance Circle: 1 : R/2 = 58
//  d = circle_inout(x+58*scale, y, 58*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//  d = circle_inout(x-58*scale, y, 58*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//
//  // Constant Reactance Circle: 1/2j : R*2 = 232
//  if (circle_inout(x, y+232*scale, 232*scale) == 0)
//    return 1;
//  if (circle_inout(x, y-232*scale, 232*scale) == 0)
//    return 1;
//
//#if 0
//  // Constant Resistance Circle: 1/3 : R*3/4 = 87
//  d = circle_inout(x+87*scale, y, 87*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//  d = circle_inout(x+87*scale, y, 87*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//#endif
//
//  // Constant Resistance Circle: 0 : R
//  d = circle_inout(x+P_RADIUS*scale, y, P_RADIUS*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//  d = circle_inout(x-P_RADIUS*scale, y, P_RADIUS*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//
//  // Constant Resistance Circle: -1/3 : R*3/2 = 174
//  d = circle_inout(x+174*scale, y, 174*scale);
//  if (d > 0) return 0;
//  if (d == 0) return 1;
//  d = circle_inout(x-174*scale, y, 174*scale);
//  //if (d > 0) return 0;
//  if (d == 0) return 1;
//  return 0;
//}

static const int cirs[][4] = {
  { 0, 58/2, 58/2, 0 },    // Constant Reactance Circle: 2j : R/2 = 58
  { 29/2, 0, 29/2, 1 },    // Constant Resistance Circle: 3 : R/4 = 29
  { 0, 115/2, 115/2, 0 },  // Constant Reactance Circle: 1j : R = 115
  { 58/2, 0, 58/2, 1 },    // Constant Resistance Circle: 1 : R/2 = 58
  { 0, 230/2, 230/2, 0 },  // Constant Reactance Circle: 1/2j : R*2 = 230
  { 86/2, 0, 86/2, 1 },    // Constant Resistance Circle: 1/3 : R*3/4 = 86
  { 0, 460/2, 460/2, 0 },  // Constant Reactance Circle: 1/4j : R*4 = 460
  { 115/2, 0, 115/2, 1 },  // Constant Resistance Circle: 0 : R
  { 173/2, 0, 173/2, 1 },  // Constant Resistance Circle: -1/3 : R*3/2 = 173
  { 0, 0, 0, 0 } // sentinel
};  

static int smith_grid3(int x, int y)
{
  int d;

  // offset to center
  x -= P_CENTER_X;
  y -= P_CENTER_Y;
  
  // outer circle
  d = circle_inout(x, y, P_RADIUS);
  if (d < 0)
    return 0;
  if (d == 0)
    return 1;

  // shift circle center to right origin
  x -= P_RADIUS /2;

  int i;
  for (i = 0; cirs[i][2]; i++) {
    d = circle_inout(x+cirs[i][0], y+cirs[i][1], cirs[i][2]);
    if (d == 0)
      return 1;
    if (d > 0 && cirs[i][3])
      return 0;
    d = circle_inout(x-cirs[i][0], y-cirs[i][1], cirs[i][2]);
    if (d == 0)
      return 1;
    if (d > 0 && cirs[i][3])
      return 0;
  }
  return 0;
}

#if 0
int
rectangular_grid(int x, int y)
{
  int c = config.grid_color;
  //#define FREQ(x) (((x) * (fspan / 1000) / (WIDTH-1)) * 1000 + fstart)
  //int32_t n = FREQ(x-1) / fgrid;
  //int32_t m = FREQ(x) / fgrid;
  //if ((m - n) > 0)
  //if (((x * 6) % (WIDTH-1)) < 6)
  //if (((x - grid_offset) % grid_width) == 0)
  if (x == 0 || x == WIDTH-1)
    return 1;
  if ((y % GRIDY) == 0)
    return 1;
  if ((((x + grid_offset) * 10) % grid_width) < 10)
    return 1;
  return 0;
}
#endif

static int rectangular_grid_x(int x)
{
  if (x < 0)
    return 0;
  if (x == 0 || x == WIDTH-1)
    return 1;
  if ((((x + grid_offset) * 10) % grid_width) < 10)
    return 1;
  return 0;
}

static int rectangular_grid_y(int y)
{
  if (y < 0)
    return 0;
  if ((y % GRIDY) == 0)
    return 1;
  return 0;
}

#if 0
int
set_strut_grid(int x)
{
  uint16_t *buf = spi_buffer;
  int y;

  for (y = 0; y < HEIGHT; y++) {
    int c = rectangular_grid(x, y);
    c |= smith_grid(x, y);
    *buf++ = 1;
  }
  return y;
}

void
draw_on_strut(int v0, int d, int color)
{
  int v;
  int v1 = v0 + d;
  if (v0 < 0) v0 = 0;
  if (v1 < 0) v1 = 0;
  if (v0 >= HEIGHT) v0 = HEIGHT-1;
  if (v1 >= HEIGHT) v1 = HEIGHT-1;
  if (v0 == v1) {
    v = v0; d = 2;
  } else if (v0 < v1) {
    v = v0; d = v1 - v0 + 1;
  } else {
    v = v1; d = v0 - v1 + 1;
  }
  while (d-- > 0)
    spi_buffer[v++] |= color;
}
#endif

/*
 * calculate log10(abs(gamma))
 */ 
static float logmag(float *v)
{
  return log10f(v[0]*v[0] + v[1]*v[1]) * 10;
}

/*
 * calculate phase[-2:2] of coefficient
 */ 
static float phase(float *v)
{
  return 2 * atan2f(v[1], v[0]) / M_PI * 90;
}

/*
 * calculate group_delay = -deltaAngle(gamma) / (deltaf * 360)
 */ 
static float group_delay(float gamma[POINT_COUNT][2], uint32_t* freq, int count, int index)
{
    float *v, *w;
    float deltaf;
    if (index == count-1) {
        deltaf = freq[index] - freq[index-1];
        v = gamma[index-1];
        w = gamma[index];
    }
    else {
        deltaf = freq[index+1] - freq[index];
        v = gamma[index];
        w = gamma[index+1];
    }
    // w = w[0]/w[1]
    // v = v[0]/v[1]
    // atan(w)-atan(v) = atan((w-v)/(1+wv))
    float r = w[0]*v[1] - w[1]*v[0];
    float i = w[0]*v[0] + w[1]*v[1];
    return atan2f(r, i) / (2 * M_PI * deltaf);
}

/*
 * calculate abs(gamma)
 */
static float linear(float *v)
{
  return  sqrtf(v[0]*v[0] + v[1]*v[1]);
}

/*
 * calculate vswr; (1+gamma)/(1-gamma)
 */ 
static float swr(float *v)
{
  float x = sqrtf(v[0]*v[0] + v[1]*v[1]);
  if (x >= 1)
    return INFINITY;
  return (1 + x)/(1 - x);
}

static float resitance(float *v) {
  float z0 = 50;
  float d = z0 / ((1-v[0])*(1-v[0])+v[1]*v[1]);
  float zr = ((1+v[0])*(1-v[0]) - v[1]*v[1]) * d;
  return zr;
}

static float reactance(float *v) {
  float z0 = 50;
  float d = z0 / ((1-v[0])*(1-v[0])+v[1]*v[1]);
  float zi = 2*v[1] * d;
  return zi;
}

static void cartesian_scale(float re, float im, int *xp, int *yp, float scale)
{
  //float scale = 4e-3;
  int x = floatToInt(re * P_RADIUS * scale);
  int y = floatToInt(im * P_RADIUS * scale);
  if (x < -P_RADIUS) x = -P_RADIUS;
  if (y < -P_RADIUS) y = -P_RADIUS;
  if (x > P_RADIUS) x = P_RADIUS;
  if (y > P_RADIUS) y = P_RADIUS;
  *xp = P_CENTER_X + x;
  *yp = P_CENTER_Y - y;
}


static uint32_t trace_into_index(
    int x, int t, int i, 
    float coeff[POINT_COUNT][2], 
    uint32_t freq[POINT_COUNT], 
    int point_count)
{
  int y = 0;
  float v = 0;
  float refpos = 10 - get_trace_refpos(t);
  float scale  = 1 / get_trace_scale(t);
  switch (trace[t].type) {
  case TRC_LOGMAG:
    v = refpos - logmag(coeff[i]) * scale;
    break;
  case TRC_PHASE:
    v = refpos - phase(coeff[i]) * scale;
    break;
  case TRC_DELAY:
    v = refpos - group_delay(coeff, freq, point_count, i) * scale;
    break;
  case TRC_LINEAR:
    v = refpos - linear(coeff[i]) * scale;
    break;
  case TRC_SWR:
    v = refpos+ (1 - swr(coeff[i])) * scale;
    break;
  case TRC_REAL:
    v = refpos - coeff[i][0] * scale;
    break;
  case TRC_IMAG:
    v = refpos - coeff[i][1] * scale;
    break;
  case TRC_R:
    v = refpos - resitance(coeff[i]) * scale;
    break;
  case TRC_X:
    v = refpos - reactance(coeff[i]) * scale;
    break;
  case TRC_SMITH:
  //case TRC_ADMIT:
  case TRC_POLAR:
    cartesian_scale(coeff[i][0], coeff[i][1], &x, &y, scale);
    return INDEX(x +CELLOFFSETX, y, i);
//    break;
  }
  if (v <  0) v =  0;
  if (v > 10) v = 10;
  y = floatToInt(v * GRIDY);
  return INDEX(x +CELLOFFSETX, y, i);
}

static int string_value_with_prefix(char *buf, int len, float val, char unit)
{
  char prefix;
  int n = 0;
  if (val < 0) {
    val = -val;
    *buf = '-';
    n++;
    len--;
  }
  if (val == INFINITY){
	  prefix = S_INFINITY[0];
  }
  else {
	  if (val < 1e-12) {
		  prefix = 'f';
		  val *= 1e15;
	  } else if (val < 1e-9) {
	  	 prefix = 'p';
	  	 val *= 1e12;
	  } else if (val < 1e-6) {
		  prefix = 'n';
		  val *= 1e9;
	  } else if (val < 1e-3) {
		  prefix = S_MICRO[0];
		  val *= 1e6;
	  } else if (val < 1) {
		  prefix = 'm';
		  val *= 1e3;
	  } else if (val < 1e3) {
		  prefix = 0;
	  } else if (val < 1e6) {
		  prefix = 'k';
		  val /= 1e3;
	  } else if (val < 1e9) {
		  prefix = 'M';
		  val /= 1e6;
	  } else {
		  prefix = 'G';
		  val /= 1e9;
	  }

	  if (val < 10) {
		  n += chsnprintf(&buf[n], len, "%.2f", val);
	  } else if (val < 100) {
		  n += chsnprintf(&buf[n], len, "%.1f", val);
	  } else {
		  n += chsnprintf(&buf[n], len, "%d", (int)val);
	  }
  }
  if (prefix)
    buf[n++] = prefix;
  if (unit)
    buf[n++] = unit;
  buf[n] = '\0';
  return n;
}


#define PI2 6.283184

static void gamma2imp(char *buf, int len, const float coeff[2], uint32_t frequency)
{
  // z = (gamma+1)/(gamma-1) * z0
  float z0 = 50;
  float d = z0 / ((1-coeff[0])*(1-coeff[0])+coeff[1]*coeff[1]);
  float zr = ((1+coeff[0])*(1-coeff[0]) - coeff[1]*coeff[1]) * d;
  float zi = 2*coeff[1] * d;
  int n;

  n = string_value_with_prefix(buf, len, zr, S_OHM[0]);
  buf[n++] = ' ';

  if (zi < 0) {
    float c = -1 / (PI2 * frequency * zi);
    string_value_with_prefix(buf+n, len-n, c, 'F');
  } else {
    float l = zi / (PI2 * frequency);
    string_value_with_prefix(buf+n, len-n, l, 'H');
  }
}

static void gamma2resistance(char *buf, int len, const float coeff[2])
{
  float z0 = 50;
  float d = z0 / ((1-coeff[0])*(1-coeff[0])+coeff[1]*coeff[1]);
  float zr = ((1+coeff[0])*(1-coeff[0]) - coeff[1]*coeff[1]) * d;
  string_value_with_prefix(buf, len, zr, S_OHM[0]);
}

static void gamma2reactance(char *buf, int len, const float coeff[2])
{
  float z0 = 50;
  float d = z0 / ((1-coeff[0])*(1-coeff[0])+coeff[1]*coeff[1]);
  float zi = 2*coeff[1] * d;
  string_value_with_prefix(buf, len, zi, S_OHM[0]);
}

static void trace_get_value_string(
    int t, char *buf, int len,
    int i, float coeff[POINT_COUNT][2], 
    uint32_t freq[POINT_COUNT], 
    int point_count)
{
  float v;
  switch (trace[t].type) {
  case TRC_LOGMAG:
    v = logmag(coeff[i]);
    if (v == -INFINITY)
      chsnprintf(buf, len, "-"S_INFINITY" dB");
    else
      chsnprintf(buf, len, "%.2fdB", v);
    break;
  case TRC_PHASE:
    v = phase(coeff[i]);
    chsnprintf(buf, len, "%.3f" S_DEGREE, v);
    break;
  case TRC_DELAY:
    v = group_delay(coeff, freq, point_count, i);
    string_value_with_prefix(buf, len, v, 's');
    break;
  case TRC_LINEAR:
    v = linear(coeff[i]);
    chsnprintf(buf, len, "%.3f", v);
    break;
  case TRC_SWR:
	v = swr(coeff[i]);
    if (v == INFINITY)
	  chsnprintf(buf, len, S_INFINITY);
    else
	  chsnprintf(buf, len, "%.3f", v);
    break;
  case TRC_SMITH:
    gamma2imp(buf, len, coeff[i], freq[i]);
    break;
  case TRC_REAL:
    chsnprintf(buf, len, "%.3f", coeff[i][0]);
    break;
  case TRC_IMAG:
    chsnprintf(buf, len, "%.3fj", coeff[i][1]);
    break;
  case TRC_R:
    gamma2resistance(buf, len, coeff[i]);
    break;
  case TRC_X:
    gamma2reactance(buf, len, coeff[i]);
    break;
  //case TRC_ADMIT:
  case TRC_POLAR:
    chsnprintf(buf, len, "%.3f %.3fj", coeff[i][0], coeff[i][1]);
    break;
  }
}

void trace_get_info(int t, char *buf, int len)
{
  const char *type = get_trace_typename(t);
  int n;
  switch (trace[t].type) {
  case TRC_LOGMAG:
    chsnprintf(buf, len, "%s %ddB/", type, (int)get_trace_scale(t));
    break;
  case TRC_PHASE:
    chsnprintf(buf, len, "%s %d" S_DEGREE "/", type, (int)get_trace_scale(t));
    break;
  case TRC_SMITH:
  //case TRC_ADMIT:
  case TRC_POLAR:
    chsnprintf(buf, len, "%s %.1fFS/", type, get_trace_scale(t));
    break;
  default:
    n = chsnprintf(buf, len, "%s ", type);
    string_value_with_prefix(buf+n, len-n, get_trace_scale(t), '/');
    break;
  }
}

static float time_of_index(int idx) {
   return 1.0 / (float)(frequencies[1] - frequencies[0]) / (float)FFT_SIZE * idx;
}

static float distance_of_index(int idx) {
#define SPEED_OF_LIGHT 299792458
   float distance = ((float)idx * (float)SPEED_OF_LIGHT) / ( (float)(frequencies[1] - frequencies[0]) * (float)FFT_SIZE * 2.0);
   return distance * (velocity_factor / 100.0);
}


static inline void mark_map(int x, int y)
{
  if (y >= 0 && y < 8 && x >= 0 && x < 16)
    markmap[current_mappage][y] |= 1<<x;
}

static inline int is_mapmarked(int x, int y)
{
  uint16_t bit = 1<<x;
  return (markmap[0][y] & bit) || (markmap[1][y] & bit);
}

static inline void markmap_upperarea(void)
{
  markmap[current_mappage][0] |= 0xffff;
}

static inline void swap_markmap(void)
{
  current_mappage = 1 - current_mappage;
}

static inline void clear_markmap(void)
{
  memset(markmap[current_mappage], 0, sizeof markmap[current_mappage]);
}

inline void force_set_markmap(void)
{
  memset(markmap[current_mappage], 0xff, sizeof markmap[current_mappage]);
}

static void mark_cells_from_index(void)
{
  int t;
  /* mark cells between each neighber points */
  for (t = 0; t < TRACE_COUNT; t++) {
    if (!trace[t].enabled)
      continue;
    int x0 = CELL_X(trace_index[t][0]);
    int y0 = CELL_Y(trace_index[t][0]);
    int m0 = x0 >> 5;
    int n0 = y0 >> 5;
    int i;
    mark_map(m0, n0);
    for (i = 1; i < sweep_points; i++) {
      int x1 = CELL_X(trace_index[t][i]);
      int y1 = CELL_Y(trace_index[t][i]);
      int m1 = x1 >> 5;
      int n1 = y1 >> 5;
      while (m0 != m1 || n0 != n1) {
        if (m0 == m1) {
          if (n0 < n1) n0++; else n0--;
        } else if (n0 == n1) {
          if (m0 < m1) m0++; else m0--;
        } else {
          int x = (m0 < m1) ? (m0 + 1)<<5 : m0<<5;
          int y = (n0 < n1) ? (n0 + 1)<<5 : n0<<5;
          int sgn = (n0 < n1) ? 1 : -1;
          if (sgn*(y-y0)*(x1-x0) < sgn*(x-x0)*(y1-y0)) {
            if (m0 < m1) m0++;
            else m0--;
          } else {
            if (n0 < n1) n0++;
            else n0--;
          }
        }
        mark_map(m0, n0);
      }
      x0 = x1;
      y0 = y1;
      m0 = m1;
      n0 = n1;
    }
  }
}

void plot_into_index(float measured[2][POINT_COUNT][2])
{
  int i, t;
  for (i = 0; i < sweep_points; i++) {
    int x = (i * (WIDTH-1) + sweep_points/2) / (sweep_points-1);
    for (t = 0; t < TRACE_COUNT; t++) {
      if (!trace[t].enabled)
        continue;
      int n = trace[t].channel;
      trace_index[t][i] = trace_into_index(
        x, t, i,
        measured[n], frequencies, sweep_points);
    }
  }
#if 0
  for (t = 0; t < TRACE_COUNT; t++)
    if (trace[t].enabled && trace[t].polar)
      quicksort(trace_index[t], 0, sweep_points);
#endif

  mark_cells_from_index();
  markmap_all_markers();
}

//static const uint8_t INSIDE = 0x00;
static const uint8_t LEFT   = 0x01;
static const uint8_t RIGHT  = 0x02;
static const uint8_t BOTTOM = 0x04;
static const uint8_t TOP    = 0x08;

static inline uint8_t _compute_outcode(int w, int h, int x, int y)
{
    uint8_t code = 0;
    if (x < 0) {
        code |= LEFT;
    } else
    if (x > w) {
        code |= RIGHT;
    }
    if (y < 0) {
        code |= BOTTOM;
    } else
    if (y > h) {
        code |= TOP;
    }
    return code;
}

static void cell_drawline(int w, int h, int x0, int y0, int x1, int y1, int c)
{
    uint8_t outcode0 = _compute_outcode(w, h, x0, y0);
    uint8_t outcode1 = _compute_outcode(w, h, x1, y1);

    if (outcode0 & outcode1) {
        // this line is out of requested area. early return
        return;
    }

    if (x0 > x1) {
        SWAP(x0, x1);
        SWAP(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    int sy = dy > 0 ? 1 : -1;
    int e = 0;

    dy *= sy;

    if (dx >= dy) {
        e = dy * 2 - dx;
        while (x0 != x1) {
            if (y0 >= 0 && y0 < h && x0 >= 0 && x0 < w)  spi_buffer[y0*w+x0] |= c;
            x0++;
            e += dy * 2;
            if (e >= 0) {
                e -= dx * 2;
                y0 += sy;
            }
        }
        if (y0 >= 0 && y0 < h && x0 >= 0 && x0 < w)  spi_buffer[y0*w+x0] |= c;
    } else {
        e = dx * 2 - dy;
        while (y0 != y1) {
            if (y0 >= 0 && y0 < h && x0 >= 0 && x0 < w)  spi_buffer[y0*w+x0] |= c;
            y0 += sy;
            e += dx * 2;
            if (e >= 0) {
                e -= dy * 2;
                x0++;
            }
        }
        if (y0 >= 0 && y0 < h && x0 >= 0 && x0 < w)  spi_buffer[y0*w+x0] |= c;
    }
}
 #if 0
int
search_index_range(int x, int y, uint32_t index[101], int *i0, int *i1)
{
  int i, j;
  int head = 0;
  int tail = sweep_points;
  i = 0;
  x &= 0x03e0;
  y &= 0x03e0;
  while (head < tail) {
    i = (head + tail) / 2;
    if (x < CELL_X0(index[i]))
      tail = i+1;
    else if (x > CELL_X0(index[i]))
      head = i;
    else if (y < CELL_Y0(index[i]))
      tail = i+1;
    else if (y > CELL_Y0(index[i]))
      head = i;
    else
      break;
  }

  if (x != CELL_X0(index[i]) || y != CELL_Y0(index[i]))
    return FALSE;
    
  j = i;
  while (j > 0 && x == CELL_X0(index[j-1]) && y == CELL_Y0(index[j-1]))
    j--;
  *i0 = j;
  j = i;
  while (j < 100 && x == CELL_X0(index[j+1]) && y == CELL_Y0(index[j+1]))
    j++;
  *i1 = j;
  return TRUE;
}
 
 #endif
 

static int search_index_range_x(int x, uint32_t index[POINT_COUNT], int *i0, int *i1)
{
  int i, j;
  int head = 0;
  int tail = sweep_points;
  x &= 0x03e0;
  i = 0;
  while (head < tail) {
    i = (head + tail) / 2;
    if (x == CELL_X0(index[i]))
      break;
    else if (x < CELL_X0(index[i])) {
      if (tail == i+1)
        break;
      tail = i+1;      
    } else {
      if (head == i)
        break;
      head = i;
    }
  }

  if (x != CELL_X0(index[i]))
    return FALSE;

  j = i;
  while (j > 0 && x == CELL_X0(index[j-1]))
    j--;
  *i0 = j;
  j = i;
  while (j < 100 && x == CELL_X0(index[j+1]))
    j++;
  *i1 = j;
  return TRUE;
}


#define X_REFERENCE_OFFSET 5
#define Y_REFERENCE_OFFSET 2
// Reference
static const uint8_t reference_bitmap[]={
0b11000000,
0b11110000,
0b11111100,
0b11110000,
0b11000000,
};
static inline void draw_refpos(int w, int h, int x, int y, int c)
{
	int y0=y;
	for (int j=0;j<5;j++,y0++)
	{
		if (y0 < 0 || y0 >= h)
			continue;
		int x0=x;
		uint8_t bits = reference_bitmap[j];
		while (bits){
			if (x0 >= 0 && x0 < w)
				spi_buffer[y0*w+x0] = (bits&0x80) ? c : DEFAULT_BG_COLOR;
			x0++;
			bits<<=1;
		}
	}
}

static void cell_draw_refpos(int m, int n, int w, int h)
{
  for (int t = 0; t < TRACE_COUNT; t++) {
    if (!trace[t].enabled || trace[t].type == TRC_SMITH || trace[t].type == TRC_POLAR)
      continue;
    int x = 0 - m * CELLWIDTH + CELLOFFSETX - X_REFERENCE_OFFSET;
    int y = 10*GRIDY - floatToInt(get_trace_refpos(t) * GRIDY) - n * CELLHEIGHT - Y_REFERENCE_OFFSET;
    if (x >=0 && x < w + 5 && y >= 0 && y < h + 5)
      draw_refpos(w, h, x, y, config.trace_color[t]);
  }
}

#define MARKER_WIDTH  7
#define MARKER_HEIGHT 10
#define X_MARKER_OFFSET 3
#define Y_MARKER_OFFSET 10
static const uint8_t marker_bitmap[]={
		// Marker 1
		0b11111110,
		0b11101110,
		0b11001110,
		0b11101110,
		0b11101110,
		0b11101110,
		0b11000110,
		0b01111100,
		0b00111000,
		0b00010000,
		// Marker 2
		0b11111110,
		0b11000110,
		0b10111010,
		0b11111010,
		0b11000110,
		0b10111110,
		0b10000010,
		0b01111100,
		0b00111000,
		0b00010000,
		// Marker 3
		0b11111110,
		0b11000110,
		0b10111010,
		0b11100110,
		0b11111010,
		0b10111010,
		0b11000110,
		0b01111100,
		0b00111000,
		0b00010000,
		// Marker 4
		0b11111110,
		0b11110110,
		0b11100110,
		0b11010110,
		0b10110110,
		0b10110110,
		0b10000010,
		0b01110100,
		0b00111000,
		0b00010000,
};

static void draw_marker(int w, int h, int x, int y, int c, int ch)
{
	int y0=y;
	for (int j=0;j<MARKER_HEIGHT;j++,y0++)
	{
		int x0=x;
		uint8_t bits = marker_bitmap[ch*10+j];
		bool force_color = false;
		while (bits){
			if (bits&0x80)
				force_color = true;
			if (x0 >= 0 && x0 < w && y0 >= 0 && y0 < h)
			{
				if (bits&0x80)
					spi_buffer[y0*w+x0] = c;
				else if (force_color)
					spi_buffer[y0*w+x0] = DEFAULT_BG_COLOR;
			}
			x0++;
			bits<<=1;
		}
	}
}

void marker_position(int m, int t, int *x, int *y)
{
    uint32_t index = trace_index[t][markers[m].index];
    *x = CELL_X(index);
    *y = CELL_Y(index);
}

int search_nearest_index(int x, int y, int t)
{
  uint32_t *index = trace_index[t];
  int min_i = -1;
  int min_d = 1000;
  int i;
  for (i = 0; i < sweep_points; i++) {
    int16_t dx = x - CELL_X(index[i]) - OFFSETX;
    int16_t dy = y - CELL_Y(index[i]) - OFFSETY;
    if (dx < 0) dx = -dx;
    if (dy < 0) dy = -dy;
    if (dx > 20 || dy > 20)
      continue;
    int d = dx*dx + dy*dy;
    if (d < min_d) {
      min_i = i;
    }
  }

  return min_i;
}

static void cell_draw_markers(int m, int n, int w, int h)
{
  int t, i;
  for (i = 0; i < MARKER_COUNT; i++) {
    if (!markers[i].enabled)
      continue;
    for (t = 0; t < TRACE_COUNT; t++) {
      if (!trace[t].enabled)
        continue;
      uint32_t index = trace_index[t][markers[i].index];
      int x = CELL_X(index) - m * CELLWIDTH - X_MARKER_OFFSET;
      int y = CELL_Y(index) - n * CELLHEIGHT - Y_MARKER_OFFSET;

      if (x >=-MARKER_WIDTH && x < w + MARKER_WIDTH && y >= -MARKER_HEIGHT && y < h + MARKER_HEIGHT)
        draw_marker(w, h, x, y, config.trace_color[t], i);
    }
  }
}

static void markmap_marker(int marker)
{
    int t;
    if (!markers[marker].enabled)
        return;
    for (t = 0; t < TRACE_COUNT; t++) {
        if (!trace[t].enabled)
            continue;
        uint32_t index = markers[marker].index;
        if (index >= POINT_COUNT)
            continue;
        index = trace_index[t][index];
        int x = CELL_X(index);
        int y = CELL_Y(index);
        int m = x>>5;
        int n = y>>5;
        mark_map(m, n);
        if ((x&31) < 6)
            mark_map(m-1, n);
        if ((x&31) > 32-6)
            mark_map(m+1, n);
        if ((y&31) < 12) {
            mark_map(m, n-1);
            if ((x&31) < 6)
                mark_map(m-1, n-1);
            if ((x&31) > 32-6)
                mark_map(m+1, n-1);
        }
    }
}

static void markmap_all_markers(void)
{
  int i;
  for (i = 0; i < MARKER_COUNT; i++) {
    if (!markers[i].enabled)
      continue;
    markmap_marker(i);
  }
  markmap_upperarea();
}


static void draw_cell(int m, int n)
{
  int x0 = m * CELLWIDTH;
  int y0 = n * CELLHEIGHT;
  int x0off = x0 - CELLOFFSETX;
  int w = CELLWIDTH;
  int h = CELLHEIGHT;
  int x, y;
  int i0, i1;
  int i;
  int t;

  if (x0off + w > area_width)
    w = area_width - x0off;
  if (y0 + h > area_height)
    h = area_height - y0;
  if (w <= 0 || h <= 0)
    return;

  chMtxLock(&mutex_ili9341); // [protect spi_buffer]
  uint16_t grid_mode = 0;
  for (t = 0; t < TRACE_COUNT; t++) {
    if (!trace[t].enabled)
      continue;

    if (trace[t].type == TRC_SMITH)
      grid_mode |= GRID_SMITH;
    //else if (trace[t].type == TRC_ADMIT)
    //  grid_mode |= GRID_ADMIT;
    else if (trace[t].type == TRC_POLAR)
      grid_mode |= GRID_POLAR;
    else
      grid_mode |= GRID_RECTANGULAR;
  }

  PULSE;
  /* draw grid */
  int c = config.grid_color;
  memset(spi_buffer, DEFAULT_BG_COLOR, sizeof spi_buffer);
  if (grid_mode & GRID_RECTANGULAR) {
    for (x = 0; x < w; x++) {
      if (rectangular_grid_x(x+x0off))
    	  for (y = 0; y < h; y++)
    		  spi_buffer[y * w + x] = c;
    }
    for (y = 0; y < h; y++) {
      if (rectangular_grid_y(y+y0))
    	  for (x = 0; x < w; x++)
    		  if (x+x0off >= 0 && x+x0off < WIDTH)
    			  spi_buffer[y * w + x] = c;
    }
  }
  if (grid_mode & (GRID_SMITH|GRID_ADMIT|GRID_POLAR)) {
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
        int n = 0;
        if (grid_mode & GRID_SMITH)
          n = smith_grid(x+x0off, y+y0);
        else if (grid_mode & GRID_ADMIT)
          n = smith_grid3(x+x0off, y+y0);
        //n = smith_grid2(x+x0, y+y0, 0.5);
        else if (grid_mode & GRID_POLAR)
          n = polar_grid(x+x0off, y+y0);
        if (n)
        	spi_buffer[y * w + x] = c;
      }
    }
  }
  PULSE;

#if 1
  /* draw rectanglar plot */
  for (t = 0; t < TRACE_COUNT; t++) {
    if (!trace[t].enabled)
      continue;
    if (trace[t].type == TRC_SMITH || trace[t].type == TRC_POLAR)
      continue;
    
    if (search_index_range_x(x0, trace_index[t], &i0, &i1)) {
      if (i0 > 0)
        i0--;
      if (i1 < POINT_COUNT-1)
        i1++;
      for (i = i0; i < i1; i++) {
        int x1 = CELL_X(trace_index[t][i]);
        int x2 = CELL_X(trace_index[t][i+1]);
        int y1 = CELL_Y(trace_index[t][i]);
        int y2 = CELL_Y(trace_index[t][i+1]);
        int c = config.trace_color[t];
        cell_drawline(w, h, x1 - x0, y1 - y0, x2 - x0, y2 - y0, c);
      }
    }
  }
#endif
#if 1
  /* draw polar plot */
  for (t = 0; t < TRACE_COUNT; t++) {
    int c = config.trace_color[t];
    if (!trace[t].enabled)
      continue;
    if (trace[t].type != TRC_SMITH && trace[t].type != TRC_POLAR)
      continue;

    for (i = 1; i < sweep_points; i++) {
      //uint32_t index = trace_index[t][i];
      //uint32_t pindex = trace_index[t][i-1];
      //if (!CELL_P(index, x0, y0) && !CELL_P(pindex, x0, y0))
      //  continue;
      int x1 = CELL_X(trace_index[t][i-1]);
      int x2 = CELL_X(trace_index[t][i]);
      int y1 = CELL_Y(trace_index[t][i-1]);
      int y2 = CELL_Y(trace_index[t][i]);
      cell_drawline(w, h, x1 - x0, y1 - y0, x2 - x0, y2 - y0, c);
    }
  }
#endif

  PULSE;
  //draw marker symbols on each trace
  cell_draw_markers(m, n, w, h);
  // draw trace and marker info on the top
  cell_draw_marker_info(m, n, w, h);
  PULSE;

  if (m == 0)
    cell_draw_refpos(m, n, w, h);

  ili9341_bulk(OFFSETX + x0off, OFFSETY + y0, w, h);
  chMtxUnlock(&mutex_ili9341); // [/protect spi_buffer]
}

static void draw_all_cells(bool flush_markmap)
{
  int m, n;
  for (m = 0; m < (area_width+CELLWIDTH-1) / CELLWIDTH; m++)
    for (n = 0; n < (area_height+CELLHEIGHT-1) / CELLHEIGHT; n++) {
      if (is_mapmarked(m, n))
        draw_cell(m, n);
    }

  if (flush_markmap) {
    // keep current map for update
    swap_markmap();
    // clear map for next plotting
    clear_markmap();
  }
}

void draw_all(bool flush)
{
    if (redraw_request & REDRAW_CELLS)
        draw_all_cells(flush);
    if (redraw_request & REDRAW_FREQUENCY)
        draw_frequencies();
    if (redraw_request & REDRAW_CAL_STATUS)
        draw_cal_status();
    redraw_request = 0;
}

void redraw_marker(int marker, int update_info)
{
  // mark map on new position of marker
  markmap_marker(marker);

  // mark cells on marker info
  if (update_info)
    markmap[current_mappage][0] = 0xffff;

  draw_all_cells(TRUE);
}

void request_to_draw_cells_behind_menu(void)
{
  int n, m;
  for (m = 7; m <= 9; m++)
    for (n = 0; n < 8; n++)
      mark_map(m, n);
  redraw_request |= REDRAW_CELLS;
}

void request_to_draw_cells_behind_numeric_input(void)
{
  int n, m;
  for (m = 0; m <= 9; m++)
    for (n = 6; n < 8; n++)
      mark_map(m, n);
  redraw_request |= REDRAW_CELLS;
}

static int cell_drawchar(int w, int h, uint8_t ch, int x, int y, int invert)
{
  uint8_t bits;
  int c, r, ch_size;
  const uint8_t *char_buf = FONT_GET_DATA(ch);
  ch_size=FONT_GET_WIDTH(ch);
  if (y <= -FONT_GET_HEIGHT || y >= h || x <= -ch_size || x >= w)
    return ch_size;
  for(c = 0; c < FONT_GET_HEIGHT; c++) {
	bits = *char_buf++;
    if ((y + c) < 0 || (y + c) >= h)
      continue;
    if (invert)
      bits = ~bits;
    for (r = 0; r < ch_size; r++) {
      if ((x+r) >= 0 && (x+r) < w && (0x80 & bits))
        spi_buffer[(y+c)*w + (x+r)] = foreground_color;
      bits <<= 1;
    }
  }
  return ch_size;
}

static void cell_drawstring(int w, int h, char *str, int x, int y)
{
  while (*str) {
    x += cell_drawchar(w, h, *str, x, y, FALSE);
    str++;
  }
}

static void cell_drawstring_invert(int w, int h, char *str, int x, int y, int invert)
{
  while (*str) {
    x += cell_drawchar(w, h, *str, x, y, invert);
    str++;
  }
}

void draw_cal_status(void)
{
#if !defined(ANTENNA_ANALYZER)
  int x = 0;
  int y = 100;
  int dy = 8;
#else
  int x = 0;
  int y = 80;
  int dy = 11;
#endif
  setForegroundColor(DEFAULT_FG_COLOR);
  setBackgroundColor(DEFAULT_BG_COLOR);
  ili9341_fill(0, y, 10, 6*dy, DEFAULT_BG_COLOR);
  if (cal_status & CALSTAT_APPLY) {
    char c[3] = "C0";
    c[1] += lastsaveid;
    if (cal_status & CALSTAT_INTERPOLATED)
      c[0] = 'c';
    else if (active_props == &current_props)
      c[1] = '*';
    ili9341_drawstring(c, x, y);
    y += dy;
  }
  if (cal_status & CALSTAT_ED) {
    ili9341_drawstring("D", x, y);
    y += dy;
  }
  if (cal_status & CALSTAT_ER) {
    ili9341_drawstring("R", x, y);
    y += dy;
  }
  if (cal_status & CALSTAT_ES) {
    ili9341_drawstring("S", x, y);
    y += dy;
  }
  if (cal_status & CALSTAT_ET) {
    ili9341_drawstring("T", x, y);
    y += dy;
  }
  if (cal_status & CALSTAT_EX) {
    ili9341_drawstring("X", x, y);
    y += dy;
  }
}

static void frequency_string(char *buf, size_t len, int32_t freq, char *prefix)
{
  if (freq < 0) {
    freq = -freq;
    *buf++ = '-';
    len -= 1;
  }
  if (freq < 1000)
  {
    chsnprintf(buf, len, "%s%d Hz", prefix, (int)freq);
  }
  else if (freq < 1000000)
  {
    chsnprintf(buf, len, "%s%d.%03d kHz", prefix,
             (int)(freq / 1000),
             (int)(freq % 1000));
  }
  else
  {
#if !defined(ANTENNA_ANALYZER)
    chsnprintf(buf, len, "%s%d.%03d %03d MHz", prefix,
             (int)(freq / 1000000),
             (int)((freq / 1000) % 1000),
             (int)(freq % 1000));
#else
	chsnprintf(buf, len, "%s%d.%03dMHz", prefix,
	         (int)(freq / 1000000),
	         (int)((freq / 1000) % 1000));
#endif
  }
}

void draw_frequencies(void)
{
  char buf1[24];
  char buf2[24];
  if ((domain_mode & DOMAIN_MODE) == DOMAIN_FREQ) {
      if (frequency1 > 0) {
        frequency_string(buf1, sizeof buf1, frequency0, "START ");
        frequency_string(buf2, sizeof buf2, frequency1, "STOP ");
      } else if (frequency1 < 0) {
        frequency_string(buf1, sizeof buf1, frequency0, "CENTER ");
        frequency_string(buf2, sizeof buf2,-frequency1, "SPAN ");
      } else {
        frequency_string(buf1, sizeof buf1, frequency0, "CW ");
      }

  } else {
	  chsnprintf(buf1, sizeof buf1, "START 0s");
      chsnprintf(buf2, sizeof buf2, "%s%d ns", "STOP ", (uint16_t)(time_of_index(POINT_COUNT-1) * 1e9));
  }
  setForegroundColor(DEFAULT_FG_COLOR);
  setBackgroundColor(DEFAULT_BG_COLOR);
#if !defined(ANTENNA_ANALYZER)
  ili9341_fill(0, 232, 320, 8, DEFAULT_BG_COLOR);
  ili9341_drawstring(buf1, OFFSETX, 232);
  ili9341_drawstring(buf2, 205, 232);
#else
  ili9341_fill(0, 231, 320, 9, DEFAULT_BG_COLOR);
  ili9341_drawstring(buf1, OFFSETX, 230);
  ili9341_drawstring(buf2, 208, 230);
#endif
}

#if !defined(ANTENNA_ANALYZER)
static void cell_draw_marker_info(int m, int n, int w, int h)
{
  char buf[24];
  int t;
  if (n != 0)
    return;
  if (active_marker < 0)
    return;
  int idx = markers[active_marker].index;
  int j = 0;
  for (t = 0; t < TRACE_COUNT; t++) {
    if (!trace[t].enabled)
      continue;
    int xpos = 1 + (j%2)*150;
    int ypos = 1 + (j/2)*8;
    xpos -= m * CELLWIDTH -CELLOFFSETX;
    ypos -= n * CELLHEIGHT;
    setForegroundColor(config.trace_color[t]);
    if (t == uistat.current_trace)
    	cell_drawstring(w, h, S_SARROW, xpos, ypos);
    xpos += 5;
    chsnprintf(buf, sizeof buf, "CH%d", trace[t].channel);
    cell_drawstring(w, h, buf, xpos, ypos);
    xpos += 19;
    trace_get_info(t, buf, sizeof buf);
    cell_drawstring(w, h, buf, xpos, ypos);
    xpos += 60;
    trace_get_value_string(
        t, buf, sizeof buf,
        idx, measured[trace[t].channel], frequencies, sweep_points);
    cell_drawstring(w, h, buf, xpos, ypos);
    j++;
  }    
  j += j&1;
  
  // LEFT
  int ypos = 1 + (j/2)*8;
  ypos -= n * CELLHEIGHT;
  if (electrical_delay != 0) {
    // draw electrical delay
    int xpos = 21;
    xpos -= m * CELLWIDTH -CELLOFFSETX;
    setForegroundColor(DEFAULT_FG_COLOR);
    chsnprintf(buf, sizeof buf, "Edelay");
    cell_drawstring(w, h, buf, xpos, ypos);
    xpos += 7 * 5;
    int n = string_value_with_prefix(buf, sizeof buf, electrical_delay * 1e-12, 's');
    cell_drawstring(w, h, buf, xpos, ypos);
    xpos += n * 5 + 5;
    float light_speed_ps = 299792458e-12; //(m/ps)
    string_value_with_prefix(buf, sizeof buf, electrical_delay * light_speed_ps * velocity_factor / 100.0, 'm');
    cell_drawstring(w, h, buf, xpos, ypos);
    ypos += 8;
  }

#ifdef __DRAW_Z__  
  {
#define ZCOLOR RGBHEX(0x00ffff)
    // draw Z
    int xpos = 1 + 2 * 5;
    xpos -= m * CELLWIDTH -CELLOFFSETX;

    float re = measured[0][idx][0];
    float im = measured[0][idx][1];
    float d = 50.0 / ((1-re)*(1-re)+im*im);
    float rs = ((1+re)*(1-re) - im*im) * d;
    float xs = 2*im * d;
    setForegroundColor(ZCOLOR);
    chsnprintf(buf, sizeof buf, "Z: %.1f %s%.1fj", rs, xs >= 0 ? "+" : "", xs);
  }
#endif

  // RIGHT
  // draw marker frequency
  int xpos = 200;
  ypos = 1 + (j/2)*8;
  xpos -= m * CELLWIDTH -CELLOFFSETX;
  ypos -= n * CELLHEIGHT;
  setForegroundColor(DEFAULT_FG_COLOR);
  chsnprintf(buf, sizeof buf, "%d:", active_marker + 1);
  cell_drawstring(w, h, buf, xpos-5-3-4, ypos);
  if ((domain_mode & DOMAIN_MODE) == DOMAIN_FREQ) {
    frequency_string(buf, sizeof buf, frequencies[idx], "");
  } else {
    //chsnprintf(buf, sizeof buf, "%d ns %.1f m", (uint16_t)(time_of_index(idx) * 1e9), distance_of_index(idx));
    int n = string_value_with_prefix(buf, sizeof buf, time_of_index(idx), 's');  
    buf[n++] = ' ';
    string_value_with_prefix(&buf[n], sizeof buf-n, distance_of_index(idx), 'm');
  }
  cell_drawstring(w, h, buf, xpos, ypos);

  // draw marker delta
  if (previous_marker >= 0 && active_marker != previous_marker && markers[previous_marker].enabled) {
    int idx0 = markers[previous_marker].index;
    xpos = 200;
    xpos -= m * CELLWIDTH -CELLOFFSETX;
    ypos += 8;
    chsnprintf(buf, sizeof buf, S_DELTA"D%d:", previous_marker+1);
    cell_drawstring(w, h, buf, xpos-7-5-5-3-4, ypos);
    if ((domain_mode & DOMAIN_MODE) == DOMAIN_FREQ) {
      frequency_string(buf, sizeof buf, frequencies[idx] - frequencies[idx0], "");
    } else {
      //chsnprintf(buf, sizeof buf, "%d ns %.1f m", (uint16_t)(time_of_index(idx) * 1e9 - time_of_index(idx0) * 1e9),
      //                                            distance_of_index(idx) - distance_of_index(idx0));
      int n = string_value_with_prefix(buf, sizeof buf, time_of_index(idx) - time_of_index(idx0), 's');
      buf[n++] = ' ';
      string_value_with_prefix(&buf[n], sizeof buf - n, distance_of_index(idx) - distance_of_index(idx0), 'm');
    }
    cell_drawstring(w, h, buf, xpos, ypos);
  }
}

#else

static void cell_draw_marker_info(int m, int n, int w, int h)
{
  char buf[24];
  int t;
  if (n > 1)
    return;
  if (active_marker < 0)
    return;
  int idx = markers[active_marker].index;
  int j = 0;
  for (t = 0; t < TRACE_COUNT; t++) {
    if (!trace[t].enabled)
      continue;
    int xpos = 1 ;
    int ypos = 1 + (j)*11;
    xpos -= m * CELLWIDTH -CELLOFFSETX;
    ypos -= n * CELLHEIGHT;

    setForegroundColor(config.trace_color[t]);
    if (t == uistat.current_trace)
    	cell_drawstring(w, h, S_SARROW, xpos, ypos);
    xpos += 6;
    chsnprintf(buf, sizeof buf, "CH%d", trace[t].channel);
    cell_drawstring(w, h, buf, xpos, ypos);

    xpos += 25;
    trace_get_info(t, buf, sizeof buf);
    cell_drawstring(w, h, buf, xpos, ypos);
    xpos += 85;
    trace_get_value_string(
        t, buf, sizeof buf,
        idx, measured[trace[t].channel], frequencies, sweep_points);
    cell_drawstring(w, h, buf, xpos, ypos);
    j++;
  }

  int ypos = 1 + (j)*11;
    ypos -= n * CELLHEIGHT;
    if (electrical_delay != 0) {
      // draw electrical delay
      int xpos = 8;
      xpos -= m * CELLWIDTH -CELLOFFSETX;
     chsnprintf(buf, sizeof buf, "Edelay");
     setForegroundColor(0xffff);
     cell_drawstring(w, h, buf, xpos, ypos);
     xpos += 7 * 7;
     int n = string_value_with_prefix(buf, sizeof buf, electrical_delay * 1e-12, 's');
     cell_drawstring(w, h, buf, xpos, ypos);
     xpos += n * 7 + 5;
     float light_speed_ps = 299792458e-12; //(m/ps)
     string_value_with_prefix(buf, sizeof buf, electrical_delay * light_speed_ps * velocity_factor / 100.0, 'm');
     cell_drawstring(w, h, buf, xpos, ypos);

   }

  // draw marker frequency
  int xpos = 195+7;
  ypos = 1;
  xpos -= m * CELLWIDTH -CELLOFFSETX;
  ypos -= n * CELLHEIGHT;
  chsnprintf(buf, sizeof buf, "%d:", active_marker + 1);
  cell_drawstring(w, h, buf, xpos, ypos);
  xpos += 12;
  if ((domain_mode & DOMAIN_MODE) == DOMAIN_FREQ) {
    frequency_string(buf, sizeof buf, frequencies[idx], "");
  } else {
    //chsnprintf(buf, sizeof buf, "%d ns %.1f m", (uint16_t)(time_of_index(idx) * 1e9), distance_of_index(idx));
    int n = string_value_with_prefix(buf, sizeof buf, time_of_index(idx), 's');  
    buf[n++] = ' ';
    string_value_with_prefix(&buf[n], sizeof buf-n, distance_of_index(idx), 'm');
  }
  cell_drawstring(w, h, buf, xpos, ypos);

  // draw marker delta
  if (previous_marker >= 0 && active_marker != previous_marker && markers[previous_marker].enabled) {
    int idx0 = markers[previous_marker].index;
    xpos = 195;
    xpos -= m * CELLWIDTH -CELLOFFSETX;
    ypos += 12;
    chsnprintf(buf, sizeof buf, "D%d:", previous_marker+1);
    cell_drawstring(w, h, buf, xpos, ypos);
    xpos += 21;
	    if ((domain_mode & DOMAIN_MODE) == DOMAIN_FREQ) {
      frequency_string(buf, sizeof buf, frequencies[idx] - frequencies[idx0], "");
    } else {
      //chsnprintf(buf, sizeof buf, "%d ns %.1f m", (uint16_t)(time_of_index(idx) * 1e9 - time_of_index(idx0) * 1e9),
      //                                            distance_of_index(idx) - distance_of_index(idx0));
      int n = string_value_with_prefix(buf, sizeof buf, time_of_index(idx) - time_of_index(idx0), 's');
      buf[n++] = ' ';
      string_value_with_prefix(&buf[n], sizeof buf - n, distance_of_index(idx) - distance_of_index(idx0), 'm');
    }
    cell_drawstring(w, h, buf, xpos, ypos);
  }
  mark_map(m, n);

}
#endif

// Draw battery level
#define BATTERY_TOP_LEVEL		4100
#define BATTERY_BOTTOM_LEVEL	3100
#define BATTERY_WARNING_LEVEL	3300

void draw_battery_status(void)
{
    uint8_t string_buf[25];
	// Set battery color
	setForegroundColor(vbat < BATTERY_WARNING_LEVEL ? RGBHEX(0xff0000) : RGBHEX(0x1fe300));
    setBackgroundColor(DEFAULT_BG_COLOR);
//    chsnprintf(string_buf, sizeof string_buf, "V:%d", vbat);
//    ili9341_drawstringV(string_buf, 1, 60);

    // Prepare battery bitmap image
    // Battery top
    int x=0;
    string_buf[x++] = 0b00111100;
    string_buf[x++] = 0b00100100;
    string_buf[x++] = 0b11111111;
//    string_buf[x++] = 0b10000001;
	// Fill battery status
	for (int power=BATTERY_TOP_LEVEL; power > BATTERY_BOTTOM_LEVEL; power-=100)
		string_buf[x++] = (power > vbat) ? 0b10000001 : // Empty line
				                           0b11111111;  // Full line
	// Battery bottom
//	string_buf[x++] = 0b10000001;
	string_buf[x++] = 0b11111111;
	// Draw battery
	blit8BitWidthBitmap(0, 1, 8, x, string_buf);
}

void
request_to_redraw_grid(void)
{
  force_set_markmap();
  redraw_request |= REDRAW_CELLS;
}

void
redraw_frame(void)
{
  ili9341_fill(0, 0, 320, 240, DEFAULT_BG_COLOR);
  draw_frequencies();
  draw_cal_status();
}

void
plot_init(void)
{
  force_set_markmap();
}
