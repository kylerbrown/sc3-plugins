/*
Non-linear Dynamic Sound Generators
Lance Putnam 2004
lance@uwalumni.com
This is a set of iterative functions and differential equations that
are known to exhibit chaotic behavior. Internal calculations are
done with 64-bit words to ensure greater accuracy.
The name of the function is followed by one of N, L, or C. These
represent the interpolation method used between function iterations.
N -> None
L -> Linear
C -> Cubic
*/
SyrinxUGen : UGen {
}
// ODEs
// 'h' is integration time-step
// Lorenz Attractor
SyrinxL : SyrinxUGen {
const <equation="x' = y\n y' = -a*g^2 - b*g^2*x - g^2*x^3 - g*x^2*y + g^2*x^2 - g*x*y";
*ar { arg freq=22050, a=0, b=0, g=24000, h=0.05, xi=0.1, yi=0, mul=1.0, add=0.0;
^this.multiNew('audio', freq, a, b, g, h, xi, yi).madd(mul, add)
}
}
