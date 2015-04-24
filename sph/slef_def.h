#ifndef _SPHDEF
#define _SPHDEF
enum Status { SOLID, LIQUID,RIGID};
#define COLOR(r,g,b)	( (unsigned int(r*255.0f)<<24) | (unsigned int(g*255.0f)<<16) | (unsigned int(b*255.0f)<<8) )
#define COLORA(r,g,b,a)	( (unsigned int(r*255.0f)<<24) | (unsigned int(g*255.0f)<<16) | (unsigned int(b*255.0f)<<8) | unsigned int(a*255.0f) )
static const float THERMAL_CONDUCTIVITY_ICE = 0.5;//0.00267;//2.18;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY_WATER = 0.1;//0.00267;//0.58;// in watts perr meter kelvin
static const float THERMAL_CONDUCTIVITY = 0.00267; //IUDN
static const float THERMAL_CONDUCTIVITY_AIR = 0.0454; //Thermal conductivity of air when it at 300 kelvin
static const float R_HEATAFFECT = 1.1;
static const float T_air = 300; //in kelvin
static const float K_WATER = 1.0; //the coefficient of interfacial tension between the water particles
static const float K_ICE = 20.0; //the coefficient of interfacial tension between water-ice particles
static const float MIN_T = 253;
static const float MAX_T = 373;
static const float ICE_T = 273;

//unsigned short m_neighbor[65536];

#endif