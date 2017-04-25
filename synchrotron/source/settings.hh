#ifndef SYNCH_SETTINGS
#define SYNCH_SETTINGS

#ifdef IO_TO_SAME_DIR
#define RESOURCE_DIR "."
#define OUTPUT_DIR "."
#else
#define RESOURCE_DIR "resources"
#define OUTPUT_DIR "calc"
#endif

namespace stron {

namespace cnst {
constexpr double pi = 3.14159265358979324;
constexpr double c = 299792458.0; // m/s
constexpr double m_proton = 938.2796e6; // eV
constexpr double s_to_turn = 11245.0;
constexpr double Qx = 0.31;
} // namespace cnst


static constexpr double FRAME_X_LOW = -2*cnst::pi;
static constexpr double FRAME_X_HIGH = 4*cnst::pi;
static constexpr double FRAME_Y_LOW = -2e9;
static constexpr double FRAME_Y_HIGH = 2e9;

static constexpr const char* META_FILE = OUTPUT_DIR"/meta.json";
static constexpr const char* PATH_FILE = OUTPUT_DIR"/particles.dat";
static constexpr const char* LINE_FILE = OUTPUT_DIR"/lines.dat";
static constexpr const char* COLL_FILE = OUTPUT_DIR"/coll.dat";
static constexpr const char* STARTDIST_FILE = OUTPUT_DIR"/startdist.dat";
static constexpr const char* ENDDIST_FILE = OUTPUT_DIR"/enddist.dat";
static constexpr const char* SIXTRACK_TEST_FILE = OUTPUT_DIR"/toymodel_track.dat";
// static constexpr const char* LHC_RAMP_FILE = "resources/ramp.txt";
static constexpr const char* LHC_RAMP_FILE = RESOURCE_DIR"/LHC_ramp.dat";
//static constexpr const char* LHC_RAMP_FILE = RESOURCE_DIR"/LHC_interpolated_60s.dat";
static constexpr const char* COLL_MOTOR_FILE = RESOURCE_DIR"/motor_tcp.txt";
//static constexpr const char* COLL_MOTOR_FILE = RESOURCE_DIR"/motor_tcp_ir3_f5433b1.txt";

} // namespace stron


#endif 
