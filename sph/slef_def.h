enum Status { SOLID, LIQUID,RIGID};
#define COLOR(r,g,b)	( (unsigned int(r*255.0f)<<24) | (unsigned int(g*255.0f)<<16) | (unsigned int(b*255.0f)<<8) )
#define COLORA(r,g,b,a)	( (unsigned int(r*255.0f)<<24) | (unsigned int(g*255.0f)<<16) | (unsigned int(b*255.0f)<<8) | unsigned int(a*255.0f) )
static const float THERMAL_CONDUCTIVITY_ICE = 0.00267;//2.18;// in watts per meter kelvin
static const float THERMAL_CONDUCTIVITY_WATER = 0.00267;//0.58;// in watts perr meter kelvin
static const float THERMAL_CONDUCTIVITY = 0.00267; //IUDN
static const float R_HEATAFFECT = 1.1;