/*
 * MagneticMapPlugin.cc — Gazebo Harmonic (gz-sim8) System Plugin
 *
 * Injects realistic magnetic anomaly + 9-term Tolles-Lawson aircraft
 * interference into the Gazebo magnetometer stream for MagNav testing.
 *
 * ── Signal model ────────────────────────────────────────────────────────────
 *
 *   Bt_total = Bt_earth(WMM)                   [from Gazebo mag sensor, nT]
 *            + B_anomaly(lat, lon)              [Ottawa map lookup, nT]
 *            + B_tl(Bx, By, Bz, coef)          [9-term TL model, nT]
 *            + noise ~ N(0, sigma)              [Gaussian noise, nT]
 *
 * The modified reading is published as a gz::msgs::Magnetometer on
 *   /magnav/magnetometer        — complete modified field (Tesla)
 *   /magnav/mag_anomaly_nT      — scalar Ottawa anomaly at current GPS
 *   /magnav/tl_interference_nT  — scalar TL aircraft interference
 *
 * PX4's gz_bridge should subscribe to /magnav/magnetometer instead of
 * the default magnetometer sensor topic.
 *
 * ── Map format (ottawa_map.bin) ──────────────────────────────────────────────
 *   uint32_t  nlat, nlon
 *   float64_t lat0_rad, dlat_rad, lon0_rad, dlon_rad, alt_m
 *   float32_t data[nlat * nlon]   — row-major, NaN = outside survey
 *
 * Built by: julia gazebo/plugins/export_ottawa_map.jl
 *
 * ── SDF parameters ───────────────────────────────────────────────────────────
 *   <map_file>    path/to/ottawa_map.bin          (required)
 *   <coef_file>   path/to/custom_tl_coef.txt      (required)
 *   <noise_nT>    5.0                              (optional, default 5 nT)
 *   <enable_tl>   true                             (optional)
 *   <enable_anomaly> true                          (optional)
 *   <world>       default                          (optional)
 *   <model>       x500_magnav                      (optional)
 *
 * ── Build ────────────────────────────────────────────────────────────────────
 *   cmake -B build gazebo/plugins && cmake --build build -j$(nproc)
 */

#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Gazebo
#include <gz/plugin/Register.hh>
#include <gz/sim/System.hh>
#include <gz/transport/Node.hh>
#include <gz/math/Vector3.hh>

// gz-msgs protobuf
#include <gz/msgs/magnetometer.pb.h>
#include <gz/msgs/navsat.pb.h>
#include <gz/msgs/float.pb.h>

namespace magnav {

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

static constexpr double kTeslaToNT  = 1.0e9;   ///< 1 T = 1e9 nT
static constexpr double kNTToTesla  = 1.0e-9;
static constexpr double kDeg2Rad    = M_PI / 180.0;

// ─────────────────────────────────────────────────────────────────────────────
// Ottawa anomaly map (loaded from binary file)
// ─────────────────────────────────────────────────────────────────────────────

/// Binary map: fixed-step regular lat/lon grid of float32 anomaly values (nT).
/// Coordinates are stored in radians (as exported by Julia).
struct OttawaMap {
  uint32_t nlat{0}, nlon{0};
  double lat0_rad{0}, dlat_rad{1};
  double lon0_rad{0}, dlon_rad{1};
  double alt_m{395.0};
  std::vector<float> data;   ///< row-major [nlat * nlon], NaN = masked

  /// Load from the binary file produced by export_ottawa_map.jl.
  bool Load(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
      gzerr << "[MagneticMapPlugin] Cannot open map file: " << path << "\n";
      return false;
    }

    // Header
    f.read(reinterpret_cast<char*>(&nlat),     4);
    f.read(reinterpret_cast<char*>(&nlon),     4);
    f.read(reinterpret_cast<char*>(&lat0_rad), 8);
    f.read(reinterpret_cast<char*>(&dlat_rad), 8);
    f.read(reinterpret_cast<char*>(&lon0_rad), 8);
    f.read(reinterpret_cast<char*>(&dlon_rad), 8);
    f.read(reinterpret_cast<char*>(&alt_m),    8);

    if (!f || nlat == 0 || nlon == 0) {
      gzerr << "[MagneticMapPlugin] Map header read failed or empty grid.\n";
      return false;
    }

    // Data
    const size_t n = static_cast<size_t>(nlat) * nlon;
    data.resize(n);
    f.read(reinterpret_cast<char*>(data.data()),
           static_cast<std::streamsize>(n * sizeof(float)));

    if (!f) {
      gzerr << "[MagneticMapPlugin] Map data truncated (expected "
            << n << " float32 values).\n";
      return false;
    }

    gzmsg << "[MagneticMapPlugin] Map loaded: " << nlat << " x " << nlon
          << "  lat=[" << lat0_rad * 180.0 / M_PI << ", "
          << (lat0_rad + (nlat - 1) * dlat_rad) * 180.0 / M_PI << "]°"
          << "  lon=[" << lon0_rad * 180.0 / M_PI << ", "
          << (lon0_rad + (nlon - 1) * dlon_rad) * 180.0 / M_PI << "]°\n";
    return true;
  }

  /// Bilinear interpolation at (lat_deg, lon_deg).
  /// Returns 0.0 outside map coverage or at masked (NaN) cells.
  double Query(double lat_deg, double lon_deg) const {
    if (data.empty()) return 0.0;

    const double lat_r = lat_deg * kDeg2Rad;
    const double lon_r = lon_deg * kDeg2Rad;

    // Fractional grid indices
    const double fi = (lat_r - lat0_rad) / dlat_rad;
    const double fj = (lon_r - lon0_rad) / dlon_rad;

    if (fi < 0.0 || fi >= nlat - 1 || fj < 0.0 || fj >= nlon - 1)
      return 0.0;   // outside coverage

    const int    i0 = static_cast<int>(fi);
    const int    j0 = static_cast<int>(fj);
    const double di = fi - i0;
    const double dj = fj - j0;

    // Four corners
    auto cell = [&](int i, int j) -> double {
      float v = data[static_cast<size_t>(i) * nlon + j];
      return std::isnan(v) ? 0.0 : static_cast<double>(v);
    };

    const double v00 = cell(i0,   j0);
    const double v01 = cell(i0,   j0+1);
    const double v10 = cell(i0+1, j0);
    const double v11 = cell(i0+1, j0+1);

    return v00 * (1-di)*(1-dj)
         + v01 * (1-di)*  dj
         + v10 *   di  *(1-dj)
         + v11 *   di  *  dj;
  }

  bool InCoverage(double lat_deg, double lon_deg) const {
    const double fi = (lat_deg * kDeg2Rad - lat0_rad) / dlat_rad;
    const double fj = (lon_deg * kDeg2Rad - lon0_rad) / dlon_rad;
    return fi >= 0.0 && fi < nlat - 1 && fj >= 0.0 && fj < nlon - 1;
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// TL coefficient loader
// ─────────────────────────────────────────────────────────────────────────────

/// Parse custom_tl_coef.txt → 9-element array (nT units).
bool LoadTLCoef(const std::string &path, std::array<double, 9> &coef) {
  std::ifstream f(path);
  if (!f) {
    gzerr << "[MagneticMapPlugin] Cannot open TL coef file: " << path << "\n";
    return false;
  }

  int idx = 0;
  std::string line;
  while (std::getline(f, line) && idx < 9) {
    if (line.size() < 2 || line[0] != 'c' || line[1] != '[') continue;
    const auto eq_pos = line.find('=');
    if (eq_pos == std::string::npos) continue;
    std::istringstream ss(line.substr(eq_pos + 1));
    double val;
    if (ss >> val) coef[idx++] = val;
  }

  if (idx != 9) {
    gzerr << "[MagneticMapPlugin] Expected 9 TL coefficients, found " << idx
          << " in " << path << "\n";
    return false;
  }

  gzmsg << "[MagneticMapPlugin] TL coef loaded: ["
        << coef[0] << ", " << coef[1] << ", " << coef[2] << ", ...]\n";
  return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// 9-term TL interference
// ─────────────────────────────────────────────────────────────────────────────

/// Compute aircraft interference (nT) using the 9-term TL design matrix.
///
/// Matches custom_tl.jl exactly:
///   cols 0-2  permanent:         Bx/Bt, By/Bt, Bz/Bt
///   cols 3-5  linear induced:    Bx,    By,    Bz
///   cols 6-8  quadratic induced: Bx²/Bt, By²/Bt, Bz²/Bt
///
/// Inputs are in nT. Returns interference in nT.
inline double TLInterference(
    double Bx_nT, double By_nT, double Bz_nT,
    const std::array<double, 9> &c)
{
  const double Bt     = std::sqrt(Bx_nT*Bx_nT + By_nT*By_nT + Bz_nT*Bz_nT);
  const double Bt_inv = 1.0 / std::max(Bt, 1.0);   // 1 nT floor

  return c[0] * (Bx_nT * Bt_inv)          // permanent
       + c[1] * (By_nT * Bt_inv)
       + c[2] * (Bz_nT * Bt_inv)
       + c[3] * Bx_nT                     // linear induced
       + c[4] * By_nT
       + c[5] * Bz_nT
       + c[6] * (Bx_nT * Bx_nT * Bt_inv) // quadratic induced
       + c[7] * (By_nT * By_nT * Bt_inv)
       + c[8] * (Bz_nT * Bz_nT * Bt_inv);
}

// ─────────────────────────────────────────────────────────────────────────────
// MagneticMapPlugin — Gazebo Harmonic System plugin
// ─────────────────────────────────────────────────────────────────────────────

class MagneticMapPlugin
    : public gz::sim::System
    , public gz::sim::ISystemConfigure
{
 public:
  MagneticMapPlugin() = default;
  ~MagneticMapPlugin() override = default;

  // ── ISystemConfigure ──────────────────────────────────────────────────────

  void Configure(
      const gz::sim::Entity & /*entity*/,
      const std::shared_ptr<const sdf::Element> &sdf,
      gz::sim::EntityComponentManager & /*ecm*/,
      gz::sim::EventManager & /*eventMgr*/) override
  {
    // ── SDF parameters ───────────────────────────────────────────────────────
    const std::string map_file  = SdfParam(sdf, "map_file",  "");
    const std::string coef_file = SdfParam(sdf, "coef_file", "");
    noise_std_nT_ = std::stod(SdfParam(sdf, "noise_nT",   "5.0"));
    enable_tl_     = SdfBool(sdf,  "enable_tl",      true);
    enable_anomaly_= SdfBool(sdf,  "enable_anomaly",  true);
    const std::string world = SdfParam(sdf, "world",  "default");
    const std::string model = SdfParam(sdf, "model",  "x500_magnav");
    const std::string link  = SdfParam(sdf, "link",   "base_link");

    // ── Load Ottawa anomaly map ───────────────────────────────────────────────
    if (map_file.empty()) {
      gzerr << "[MagneticMapPlugin] <map_file> SDF parameter is required.\n";
      return;
    }
    if (!map_.Load(map_file)) return;

    // ── Load TL coefficients ─────────────────────────────────────────────────
    if (enable_tl_) {
      if (coef_file.empty()) {
        gzerr << "[MagneticMapPlugin] <coef_file> required when enable_tl=true.\n";
        return;
      }
      if (!LoadTLCoef(coef_file, coef_)) return;
    }

    gzmsg << "[MagneticMapPlugin] noise=" << noise_std_nT_ << " nT"
          << "  TL=" << (enable_tl_ ? "on" : "off")
          << "  anomaly=" << (enable_anomaly_ ? "on" : "off") << "\n";

    // ── Topic names ───────────────────────────────────────────────────────────
    // NavSat topic: constructed from world/model/link (not overridable — there
    // is only one GPS sensor and its topic follows the standard Gazebo pattern).
    navsat_topic_ = "/world/" + world + "/model/" + model +
                    "/link/" + link + "/sensor/navsat_sensor/navsat";

    // Magnetometer INPUT topic:
    //   If <mag_input_topic> is set in SDF, use that (e.g. /magnav/raw_magnetometer
    //   when the sensor's <topic> is redirected in model.sdf).
    //   Default: the standard Gazebo sensor topic (for non-redirected models).
    const std::string default_mag_in = "/world/" + world + "/model/" + model +
                    "/link/" + link + "/sensor/magnetometer_sensor/magnetometer";
    mag_in_topic_ = SdfParam(sdf, "mag_input_topic", default_mag_in);

    // Magnetometer OUTPUT topic:
    //   If <mag_output_topic> is set, publish corrected field there.
    //   Default: /magnav/magnetometer (monitoring only; PX4 won't see it unless
    //   PX4's gz_bridge is told to subscribe here).
    //   For the full no-PX4-changes setup, set this to the PX4-expected topic:
    //     /world/{world}/model/{model}/link/{link}/sensor/magnetometer_sensor/magnetometer
    mag_out_topic_ = SdfParam(sdf, "mag_output_topic", "/magnav/magnetometer");

    // ── Subscribe to NavSat (GPS) ─────────────────────────────────────────────
    if (!node_.Subscribe(navsat_topic_,
                         &MagneticMapPlugin::OnNavSat, this)) {
      gzwarn << "[MagneticMapPlugin] Could not subscribe to " << navsat_topic_
             << " — will use default Ottawa position\n";
    } else {
      gzmsg << "[MagneticMapPlugin] Subscribed GPS: " << navsat_topic_ << "\n";
    }

    // ── Subscribe to raw magnetometer ────────────────────────────────────────
    if (!node_.Subscribe(mag_in_topic_,
                         &MagneticMapPlugin::OnMagnetometer, this)) {
      gzerr << "[MagneticMapPlugin] Could not subscribe to " << mag_in_topic_ << "\n";
      return;
    }
    gzmsg << "[MagneticMapPlugin] Subscribed Mag (in):  " << mag_in_topic_  << "\n";
    gzmsg << "[MagneticMapPlugin] Publishing Mag (out): " << mag_out_topic_ << "\n";

    // ── Advertise output topics ───────────────────────────────────────────────
    // Primary output: corrected field on the configured output topic.
    // This is what PX4 gz_bridge reads when mag_output_topic is set to the
    // PX4-expected sensor topic.
    pub_mag_     = node_.Advertise<gz::msgs::Magnetometer>(mag_out_topic_);

    // Monitoring topics: always published, world/model-independent.
    pub_anomaly_ = node_.Advertise<gz::msgs::Float>("/magnav/mag_anomaly_nT");
    pub_tl_      = node_.Advertise<gz::msgs::Float>("/magnav/tl_interference_nT");

    gzmsg << "[MagneticMapPlugin] Publishing /magnav/mag_anomaly_nT\n";
    gzmsg << "[MagneticMapPlugin] Publishing /magnav/tl_interference_nT\n";
    gzmsg << "[MagneticMapPlugin] Ready.\n";
  }

 private:
  // ── GPS callback ──────────────────────────────────────────────────────────
  void OnNavSat(const gz::msgs::NavSat &msg)
  {
    std::lock_guard<std::mutex> lock(state_mutex_);
    lat_deg_   = msg.latitude_deg();
    lon_deg_   = msg.longitude_deg();
    gps_valid_ = true;
  }

  // ── Magnetometer callback — core computation ──────────────────────────────
  void OnMagnetometer(const gz::msgs::Magnetometer &msg)
  {
    // Snapshot GPS state
    double lat_deg, lon_deg;
    {
      std::lock_guard<std::mutex> lock(state_mutex_);
      lat_deg = lat_deg_;
      lon_deg = lon_deg_;
    }

    // Raw Earth field in Tesla from Gazebo (WMM model)
    const double Bx_T = msg.field_tesla().x();
    const double By_T = msg.field_tesla().y();
    const double Bz_T = msg.field_tesla().z();

    // Convert to nT for all calculations
    const double Bx_nT = Bx_T * kTeslaToNT;
    const double By_nT = By_T * kTeslaToNT;
    const double Bz_nT = Bz_T * kTeslaToNT;
    const double Bt_earth_nT = std::sqrt(Bx_nT*Bx_nT + By_nT*By_nT + Bz_nT*Bz_nT);

    // 1) Ottawa anomaly map lookup
    const double anomaly_nT = enable_anomaly_ ? map_.Query(lat_deg, lon_deg) : 0.0;

    // 2) 9-term TL aircraft interference
    const double tl_nT = enable_tl_
        ? TLInterference(Bx_nT, By_nT, Bz_nT, coef_)
        : 0.0;

    // 3) Gaussian noise  (thread_local rng avoids mutex on the hot path)
    thread_local std::mt19937 rng{std::random_device{}()};
    thread_local std::normal_distribution<double> dist;
    const double noise_nT = dist(rng, std::normal_distribution<double>::param_type{0.0, noise_std_nT_});

    // 4) Total scalar field
    const double Bt_total_nT = Bt_earth_nT + anomaly_nT + tl_nT + noise_nT;

    // Scale the direction vector so its magnitude equals Bt_total.
    // (Preserves field direction from WMM; adjusts magnitude only.)
    const double scale = (Bt_earth_nT > 1.0) ? (Bt_total_nT / Bt_earth_nT) : 1.0;

    // ── Publish /magnav/magnetometer ─────────────────────────────────────────
    gz::msgs::Magnetometer out_mag;
    *out_mag.mutable_header() = msg.header();   // copy timestamp
    out_mag.mutable_field_tesla()->set_x(Bx_T * scale);
    out_mag.mutable_field_tesla()->set_y(By_T * scale);
    out_mag.mutable_field_tesla()->set_z(Bz_T * scale);
    pub_mag_.Publish(out_mag);

    // ── Publish anomaly scalar ────────────────────────────────────────────────
    gz::msgs::Float f_anomaly;
    f_anomaly.set_data(static_cast<float>(anomaly_nT));
    pub_anomaly_.Publish(f_anomaly);

    // ── Publish TL scalar ─────────────────────────────────────────────────────
    gz::msgs::Float f_tl;
    f_tl.set_data(static_cast<float>(tl_nT));
    pub_tl_.Publish(f_tl);

    // Debug log (throttled to every 500 messages ≈ every 5 s at 100 Hz)
    if ((msg_count_++ % 500) == 0) {
      gzmsg << "[MagneticMapPlugin] "
            << "lat=" << lat_deg << "° lon=" << lon_deg << "°"
            << (map_.InCoverage(lat_deg, lon_deg) ? " [IN map]" : " [OUT of map]")
            << "  earth=" << static_cast<int>(Bt_earth_nT)
            << "  anom=" << static_cast<int>(anomaly_nT)
            << "  tl=" << static_cast<int>(tl_nT)
            << "  total=" << static_cast<int>(Bt_total_nT) << " nT"
            << "  → " << mag_out_topic_ << "\n";
    }
  }

  // ── SDF helpers ───────────────────────────────────────────────────────────
  static std::string SdfParam(
      const std::shared_ptr<const sdf::Element> &sdf,
      const std::string &key, const std::string &def)
  {
    if (!sdf || !sdf->HasElement(key)) return def;
    return sdf->Get<std::string>(key);
  }
  static bool SdfBool(
      const std::shared_ptr<const sdf::Element> &sdf,
      const std::string &key, bool def)
  {
    if (!sdf || !sdf->HasElement(key)) return def;
    return sdf->Get<bool>(key);
  }

  // ── State ─────────────────────────────────────────────────────────────────
  OttawaMap            map_;
  std::array<double,9> coef_{};
  double               noise_std_nT_{5.0};
  bool                 enable_tl_{true};
  bool                 enable_anomaly_{true};

  gz::transport::Node             node_;
  gz::transport::Node::Publisher  pub_mag_;
  gz::transport::Node::Publisher  pub_anomaly_;
  gz::transport::Node::Publisher  pub_tl_;

  std::string navsat_topic_;
  std::string mag_in_topic_;    ///< raw WMM field input (redirected sensor topic)
  std::string mag_out_topic_;   ///< corrected field output (PX4 or monitor topic)

  mutable std::mutex state_mutex_;
  double  lat_deg_{44.76};   // Ottawa default until first GPS fix
  double  lon_deg_{-76.15};
  bool    gps_valid_{false};

  std::atomic<uint64_t> msg_count_{0};
};

}  // namespace magnav

// ── Plugin registration ───────────────────────────────────────────────────────
GZ_ADD_PLUGIN(magnav::MagneticMapPlugin,
              gz::sim::System,
              gz::sim::ISystemConfigure)

GZ_ADD_PLUGIN_ALIAS(magnav::MagneticMapPlugin,
                    "magnav::MagneticMapPlugin")
