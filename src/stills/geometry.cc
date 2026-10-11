#include "geometry.hh"

#include <cmath>

namespace mxi::stills {

double StillsPanel::attenuation_length(const Vec3 &u) const {
  // The normal is f x s unnormalised, turned to face away from the sample.
  Vec3 normal = fast.cross(slow);
  if (origin.dot(normal) < 0.0)
    normal = -normal;
  const double cos_t = u.dot(normal);
  if (!(mu > 0.0 && cos_t > 0.0))
    throw GeometryError("parallax correction needs mu > 0 and a ray meeting "
                        "the sensor from the front; got mu = " +
                        std::to_string(mu) +
                        ", cos_t = " + std::to_string(cos_t));
  return (1.0 / mu) - (t0 / cos_t + 1.0 / mu) * std::exp(-mu * t0 / cos_t);
}

std::pair<double, double> StillsPanel::px_to_mm(double px_fast,
                                                double px_slow) const {
  const double x = px_fast * pixel_size[0];
  const double y = px_slow * pixel_size[1];
  if (!parallax())
    return {x, y};
  const Vec3 u = (origin + fast * x + slow * y).normalized();
  const double length = attenuation_length(u);
  return {x - u.dot(fast) * length, y - u.dot(slow) * length};
}

const StillsPanel &StillsDetector::panel(std::int64_t index) const {
  if (index < 0 || static_cast<std::size_t>(index) >= panels.size())
    throw GeometryError("panel " + std::to_string(index) +
                        " is outside a detector of " +
                        std::to_string(panels.size()) + " panels");
  return panels[static_cast<std::size_t>(index)];
}

MappedSpots map_spots(const StillsDetector &detector, const Vec3 &s0,
                      double wavelength, const Spots &spots) {
  const std::size_t n = spots.size();
  if (spots.x_px.size() != n || spots.y_px.size() != n ||
      spots.var_x_px.size() != n || spots.var_y_px.size() != n)
    throw GeometryError("spot columns differ in length");

  MappedSpots out;
  out.x_mm.resize(n);
  out.y_mm.resize(n);
  out.var_x_mm.resize(n);
  out.var_y_mm.resize(n);
  out.s1.resize(n);
  out.rlp.resize(n);
  const double inverse_wavelength = 1.0 / wavelength;
  for (std::size_t i = 0; i < n; ++i) {
    const StillsPanel &p = detector.panel(spots.panel[i]);
    const auto [x, y] = p.px_to_mm(spots.x_px[i], spots.y_px[i]);
    out.x_mm[i] = x;
    out.y_mm[i] = y;
    out.var_x_mm[i] = spots.var_x_px[i] * (p.pixel_size[0] * p.pixel_size[0]);
    out.var_y_mm[i] = spots.var_y_px[i] * (p.pixel_size[1] * p.pixel_size[1]);
    const Vec3 lab = p.lab_coord(x, y);
    const Vec3 s1 = (lab / lab.norm()) * inverse_wavelength;
    out.s1[i] = s1;
    out.rlp[i] = s1 - s0;
  }
  return out;
}

} // namespace mxi::stills
