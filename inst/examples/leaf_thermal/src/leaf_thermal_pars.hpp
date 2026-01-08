namespace odelia {
namespace ode {

struct LeafThermalPars
{
  // Air temperature sinusoid
  double Tmean = 32.0; // °C
  double Tamp = 6.0;   // °C amplitude
  double tpeak = 15.0; // hour of daily peak

  // Leaf thermal behaviour
  double k_H = 0.5;       // 1/h
  double g_tr_max = 1.0;  // °C/h cooling
  double m_tr = 0.5;      // steepness
  double T_tr_mid = 30.0; // °C
};
}
}
