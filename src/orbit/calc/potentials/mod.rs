//! This module provides realizations for
//! different models of Galactic potential
//!
//! Here's a list of available potentials and their R and Z derivatives.
//! Note that $ r^2 = X^2 + Y^2 + Z^2 = R^2 + Z^2 $.
//!
//! 1. Plummer (P) potential with parameters $ (M, b) $:
//!
//! $$
//! \Phi(r(R, Z)) = - \frac{M}{(r^2 + b^2)^{1/2}};
//! $$
//!
//! $$
//! \frac{\partial \Phi(R, Z)}{\partial R} =
//! \frac{M R}{(R^2 + Z^2 + b^2)^{3/2}};
//! $$
//!
//! $$
//! \frac{\partial \Phi(R, Z)}{\partial Z} =
//! \frac{M Z}{(R^2 + Z^2 + b^2)^{3/2}}.
//! $$
//!
//! 2. Miyamoto & Nagai (MN) potential with parameters $ (M, a, b) $:
//!
//! $$
//! \Phi(R, Z) = - \frac{M}{\left[ R^2 +
//! \left( a + \sqrt{Z^2 + b^2} \right)^2 \right]^{1/2}};
//! $$
//!
//! $$
//! \frac{\partial \Phi(R, Z)}{\partial R} =
//! \frac{M R}{\left[ R^2 + \left( a +
//! \sqrt{Z^2 + b^2} \right)^2 \right]^{3/2}};
//! $$
//!
//! $$
//! \frac{\partial \Phi(R, Z)}{\partial Z} =
//! \frac{M Z \left( a + \sqrt{b^2 + Z^2} \right)}{
//! \sqrt{b^2 + Z^2} \left( R^2 + \left( a +
//! \sqrt{b^2 + Z^2} \right)^2 \right)^{3/2}}.
//! $$
//!
//! 3. Navarro-Frenk-White (NFW) potential with parameters $ (M, a) $:
//!
//! $$
//! \Phi(r(R, Z)) = - \frac{M}{r} \ln{\left( 1 +
//! \frac{r}{a} \right)};
//! $$
//!
//! $$
//! \frac{\partial \Phi(R, Z)}{\partial R} =
//! \frac{M R \ln{\left( \sqrt{R^2 + Z^2} / a
//! + 1 \right)}}{\left( R^2 + Z^2 \right)^{3/2}} -
//! \frac{M R}{\left( R^2 + Z^2 \right) \left(
//! \sqrt{R^2+ Z^2} + a \right)};
//! $$
//!
//! $$
//! \frac{\partial \Phi(R, Z)}{\partial Z} =
//! \frac{M Z \ln{\left( \sqrt{R^2 + Z^2} / a
//! + 1 \right)}}{\left( R^2 + Z^2 \right)^{3/2}} -
//! \frac{M Z}{\left( R^2 + Z^2 \right) \left(
//! \sqrt{R^2+ Z^2} + a \right)};
//! $$

pub mod miyamoto_nagai;
pub mod navarro_frenk_white;
pub mod plummer;
