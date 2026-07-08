#ifndef APPARENT_HORIZON_FINDER_AHF_DATA_H_
#define APPARENT_HORIZON_FINDER_AHF_DATA_H_

namespace exahype2::ccz4::ahf {

/**
 * Cartesian reference coordinate for one apparent-horizon centre point.
 */
struct ReferencePoint {
  double x;
  double y;
  double z;
};

}  // namespace exahype2::ccz4::ahf

#endif  // APPARENT_HORIZON_FINDER_AHF_DATA_H_
