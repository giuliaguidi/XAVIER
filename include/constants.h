
#define NINF  	(std::numeric_limits<int8_t>::min())
#define goRIGHT (0)
#define goDOWN  (1)
#define MIDDLE 	(LOGICALWIDTH / 2)

#define CUTOFF	(std::numeric_limits<int8_t>::max() - 25)

#ifdef DEBUG
	#define myLog( var ) do { std::cerr << "LOG:	" << __FILE__ << "(" << __LINE__ << ")	" << #var << " = " << (var) << std::endl; } while(0)
#else
	#define myLog( var )
#endif


enum ExtDirectionX
{
	XAVIER_EXTEND_NONE  = 0,
	XAVIER_EXTEND_LEFT  = 1,
	XAVIER_EXTEND_RIGHT = 2,
	XAVIER_EXTEND_BOTH  = 3
};