// Define a struct to hold parameters for the leaf thermal model
struct LeafThermalPars
{
  // Leaf thermal behaviour
  double k_H = 0.5;       // 1/h
  double g_tr_max = 1.0;  // °C/h cooling
  double m_tr = 0.5;      // steepness
  double T_tr_mid = 30.0; // °C
};
