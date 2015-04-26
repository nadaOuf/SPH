#ifndef _SPHDEF
#define _SPHDEF
enum Status {SOLID, LIQUID, RIGID, SOURCE};
#define COLOR(r,g,b)	( (unsigned int(r*255.0f)<<24) | (unsigned int(g*255.0f)<<16) | (unsigned int(b*255.0f)<<8) )
#define COLORA(r,g,b,a)	( (unsigned int(r*255.0f)<<24) | (unsigned int(g*255.0f)<<16) | (unsigned int(b*255.0f)<<8) | unsigned int(a*255.0f) )

static const float THERMAL_CONDUCTIVITY_ICE = 0.00267f;//2.18;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY_WATER = 0.00267f;//0.58;// in watts perr meter kelvin
static const float THERMAL_CONDUCTIVITY = 0.00267f; //IUDN
static const float THERMAL_CONDUCTIVITY_AIR = 0.0454f; //Thermal conductivity of air when it at 300 kelvin
static const float R_HEATAFFECT = 1.1f;
static const float T_air = 300.0f; //in kelvin
static const float K_WATER = 1.0; //the coefficient of interfacial tension between the water particles
static const float K_ICE = 20.0f; //the coefficient of interfacial tension between water-ice particles
static const float MIN_T = 253.0f;
static const float MAX_T = 373.0f;
static const float ICE_T = 273.0f;
static const float HEAT_CAPACITY_ICE = 2.03f;
static const float HEAT_CAPACITY_WATER = 4.179f;
static const float Bounce =1.65;
static const float EPSILON = 0.01;

static const float Kw = 0.5f;
static const float Kice = 2.0f;

//unsigned short m_neighbor[65536];

#endif
