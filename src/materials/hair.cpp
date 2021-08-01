
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// materials/hair.cpp*
#include <array>
#include <numeric>
#include "interaction.h"
#include "materials/hair.h"
#include "paramset.h"
#include "reflection.h"
#include "sampling.h"
#include "spectrum.h"
#include "texture.h"
#include "textures/constant.h"

// (Mandy Xia) Add this to be able to acquire curve primitives info.
#include "../shapes/curve.h"
extern float *BSDFTABLE;
extern float *PDFTABLE;
extern float *CDFTABLE;
extern float *RATIO;

namespace pbrt {

// Hair Local Declarations
inline Float I0(Float x), LogI0(Float x);

// Hair Local Functions
static Float Mp(Float cosThetaI, Float cosThetaO, Float sinThetaI,
                Float sinThetaO, Float v) {
    Float a = cosThetaI * cosThetaO / v;
    Float b = sinThetaI * sinThetaO / v;
    Float mp =
        (v <= .1)
            ? (std::exp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
            : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
    CHECK(!std::isinf(mp) && !std::isnan(mp));
    return mp;
}

inline Float I0(Float x) {
    Float val = 0;
    Float x2i = 1;
    int64_t ifact = 1;
    int i4 = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for (int i = 0; i < 10; ++i) {
        if (i > 1) ifact *= i;
        val += x2i / (i4 * Sqr(ifact));
        x2i *= x * x;
        i4 *= 4;
    }
    return val;
}

inline Float LogI0(Float x) {
    if (x > 12)
        return x + 0.5 * (-std::log(2 * Pi) + std::log(1 / x) + 1 / (8 * x));
    else
        return std::log(I0(x));
}

static std::array<Spectrum, pMax + 1> Ap(Float cosThetaO, Float eta, Float h,
                                         const Spectrum &T) {
    std::array<Spectrum, pMax + 1> ap;
    // Compute $p=0$ attenuation at initial cylinder intersection
    Float cosGammaO = SafeSqrt(1 - h * h);
    Float cosTheta = cosThetaO * cosGammaO;
    Float f = FrDielectric(cosTheta, 1.f, eta);
    ap[0] = f;

    // Compute $p=1$ attenuation term
    ap[1] = Sqr(1 - f) * T;

    // Compute attenuation terms up to $p=_pMax_$
    for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

    // Compute attenuation term accounting for remaining orders of scattering
    ap[pMax] = ap[pMax - 1] * f * T / (Spectrum(1.f) - T * f);
    return ap;
}

inline Float Phi(int p, Float gammaO, Float gammaT) {
    return 2 * p * gammaT - 2 * gammaO + p * Pi;
}

inline Float Logistic(Float x, Float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
}

inline Float LogisticCDF(Float x, Float s) {
    return 1 / (1 + std::exp(-x / s));
}

inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
    CHECK_LT(a, b);
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

inline Float Np(Float phi, int p, Float s, Float gammaO, Float gammaT) {
    Float dphi = phi - Phi(p, gammaO, gammaT);
    // Remap _dphi_ to $[-\pi,\pi]$
    while (dphi > Pi) dphi -= 2 * Pi;
    while (dphi < -Pi) dphi += 2 * Pi;
    return TrimmedLogistic(dphi, s, -Pi, Pi);
}

static Float SampleTrimmedLogistic(Float u, Float s, Float a, Float b) {
    CHECK_LT(a, b);
    Float k = LogisticCDF(b, s) - LogisticCDF(a, s);
    Float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
    CHECK(!std::isnan(x));
    return Clamp(x, a, b);
}

// HairMaterial Method Definitions
void HairMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                              MemoryArena &arena,
                                              TransportMode mode,
                                              bool allowMultipleLobes) const {
    // (Mandy Xia) Read fiber orientation from curve primitives.
    const Curve* cur = dynamic_cast<const Curve*>(si->shape);
    Float ori = 0;
    if (cur){
      ori = cur->ori;
    }

    Float bm = beta_m->Evaluate(*si);
    Float bn = beta_n->Evaluate(*si);
    Float a = alpha->Evaluate(*si);
    Float e = eta->Evaluate(*si);

    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, e);

    Spectrum sig_a;
    if (sigma_a)
        sig_a = sigma_a->Evaluate(*si).Clamp();
    else if (color) {
        Spectrum c = color->Evaluate(*si).Clamp();
        sig_a = HairBSDF::SigmaAFromReflectance(c, bn);
    } else {
        CHECK(eumelanin || pheomelanin);
        sig_a = HairBSDF::SigmaAFromConcentration(
            std::max(Float(0), eumelanin ? eumelanin->Evaluate(*si) : 0),
            std::max(Float(0), pheomelanin ? pheomelanin->Evaluate(*si) : 0));
    }

    // Offset along width
    Float h = -1 + 2 * si->uv[1];
    int wavelengthindex = si->wavelengthindex;   // (Mandy Xia) Get wavelength index from SurfaceInteraction.
    si->bsdf->Add(ARENA_ALLOC(arena, HairBSDF)(h, e, sig_a, bm, bn, a, crosssection, ori, wavelengthindex));
}

HairMaterial *CreateHairMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> sigma_a =
        mp.GetSpectrumTextureOrNull("sigma_a");
    std::shared_ptr<Texture<Spectrum>> color =
        mp.GetSpectrumTextureOrNull("color");
    std::shared_ptr<Texture<Float>> eumelanin =
        mp.GetFloatTextureOrNull("eumelanin");
    std::shared_ptr<Texture<Float>> pheomelanin =
        mp.GetFloatTextureOrNull("pheomelanin");
    if (sigma_a) {
        if (color)
            Warning(
                "Ignoring \"color\" parameter since \"sigma_a\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"sigma_a\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"sigma_a\" was "
                "provided.");
    } else if (color) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since \"color\" was provided.");
        if (eumelanin)
            Warning(
                "Ignoring \"eumelanin\" parameter since \"color\" was "
                "provided.");
        if (pheomelanin)
            Warning(
                "Ignoring \"pheomelanin\" parameter since \"color\" was "
                "provided.");
    } else if (eumelanin || pheomelanin) {
        if (sigma_a)
            Warning(
                "Ignoring \"sigma_a\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
        if (color)
            Warning(
                "Ignoring \"color\" parameter since "
                "\"eumelanin\"/\"pheomelanin\" was provided.");
    } else {
        // Default: brown-ish hair.
        sigma_a = std::make_shared<ConstantTexture<Spectrum>>(
            HairBSDF::SigmaAFromConcentration(1.3, 0.));
    }

    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTexture("eta", 1.55f);
    std::shared_ptr<Texture<Float>> beta_m = mp.GetFloatTexture("beta_m", 0.3f);
    std::shared_ptr<Texture<Float>> beta_n = mp.GetFloatTexture("beta_n", 0.3f);
    std::shared_ptr<Texture<Float>> alpha = mp.GetFloatTexture("alpha", 2.f);

    // (Mandy Xia) Read in cross-section type;
    int crosssection = mp.FindInt("crosssection", 0);
    return new HairMaterial(sigma_a, color, eumelanin, pheomelanin, eta, beta_m,
                            beta_n, alpha, crosssection);
}

// HairBSDF Method Definitions
HairBSDF::HairBSDF(Float h, Float eta, const Spectrum &sigma_a, Float beta_m,
                   Float beta_n, Float alpha, int crosssection, Float ori, int wavelengthindex)
    : BxDF(BxDFType(BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION)),
      h(h),
      gammaO(SafeASin(h)),
      eta(eta),
      sigma_a(sigma_a),
      beta_m(beta_m),
      beta_n(beta_n),
      crosssection(crosssection),
      ori(ori),
      wavelengthindex(wavelengthindex){
    CHECK(h >= -1 && h <= 1);
    CHECK(beta_m >= 0 && beta_m <= 1);
    CHECK(beta_n >= 0 && beta_n <= 1);
    // Compute longitudinal variance from $\beta_m$
    static_assert(
        pMax >= 3,
        "Longitudinal variance code must be updated to handle low pMax");
    v[0] = Sqr(0.726f * beta_m + 0.812f * Sqr(beta_m) + 3.7f * Pow<20>(beta_m));
    // (Mandy Xia) No cuticle scales, v sets to the same.
    for (int p = 1; p<=pMax;++p)
      v[p] = v[0];

    // Compute azimuthal logistic scale factor from $\beta_n$
    s = SqrtPiOver8 *
        (0.265f * beta_n + 1.194f * Sqr(beta_n) + 5.372f * Pow<22>(beta_n));
    CHECK(!std::isnan(s));

    // Compute $\alpha$ terms for hair scales
    sin2kAlpha[0] = std::sin(Radians(alpha));
    cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
    for (int i = 1; i < 3; ++i) {
        sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
        cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
    }
}


  // (Mandy Xia) Evaluate f using precomputation tables.
Spectrum HairBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
    // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);
    phiO += ori * M_PI/180;

    // Compute hair coordinate system terms related to _wi_
    Float sinThetaI = wi.x;
    Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    Float phiI = std::atan2(wi.z, wi.y);
    phiI += ori * M_PI / 180.0;  // Take orientation into account.

    // resolution
    int thetanum = 100; // longitudinal angle resolution
    int phiinum = 360;  // outgoing phi resolution
    int phionum;        // incoming phi resolution

    // Evaluate hair BSDF
    Float phi = phiI - phiO;

    Float asino = SafeASin(std::abs(sinThetaO));
    Float unit = Pi / 2.0 / (Float)thetanum;
    int floortheta = Clamp(std::floor(asino / unit), 0, thetanum-1);

    int floorphio, floorphii;
    Float phiOnew = phiO;
    while (phiOnew > 2*Pi) phiOnew -= 2 * Pi;
    while (phiOnew < 0) phiOnew += 2 * Pi;
    Float phiInew = phiI;
    while (phiInew > 2*Pi) phiInew -= 2 * Pi;
    while (phiInew < 0) phiInew += 2 * Pi;
    Float phiunit = M_PI * 2 / (Float)phiinum;
    if (crosssection==0){
      // circle
      phionum = 1;
      floorphio = 0;
      while (phi > 2*Pi) phi -= 2 * Pi;
      while (phi < 0) phi += 2 * Pi;
      floorphii = Clamp(std::floor(phi / phiunit), 0, phiinum-1);
    }else if (crosssection==1){
      // ellipse
      phionum = 90;
      //if (phiOnew < Pi/2) // don't need to change
      if (phiOnew >= Pi/2 && phiOnew < Pi){
        phiOnew = Pi - phiOnew;
        phiInew = Pi - phiInew;
        while (phiInew < 0) phiInew += 2 * Pi;
      } else if (phiOnew >= Pi && phiOnew < 3.0*Pi/2.0){
        phiOnew -= Pi;
        phiInew -= Pi;
        while (phiInew < 0) phiInew += 2 * Pi;
      } else if (phiOnew >= 3.0*Pi/2.0){
        phiOnew -= Pi;
        phiInew -= Pi;
        // treat it as second quadrant
        phiOnew = Pi - phiOnew;
        phiInew = Pi - phiInew;
        while (phiInew < 0) phiInew += 2 * Pi;
      }
      Float phiounit = Pi/2.0/(Float)phionum;
      floorphio = Clamp(std::floor(phiOnew / phiounit), 0, phionum-1);
      floorphii = Clamp(std::floor(phiInew / phiunit), 0, phiinum-1);
    }else{
      // nonelliptical
      phionum = 360;
      floorphio = Clamp(std::floor(phiOnew / phiunit), 0, phionum-1);
      floorphii = Clamp(std::floor(phiInew / phiunit), 0, phionum-1);
    }

    unsigned long long num_samples = thetanum * phionum * phiinum;
    unsigned long long num_nor = phionum * phiinum;

    Float AN = BSDFTABLE[wavelengthindex*num_samples + floortheta*num_nor + floorphio*phiinum + floorphii];

    Spectrum fsum(0.);
    cosThetaO = std::abs(cosThetaO);
    fsum = Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[0]) * AN;
    if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);
    CHECK(!std::isinf(fsum.y()) && !std::isnan(fsum.y()));

    return fsum;
}

std::array<Float, pMax + 1> HairBSDF::ComputeApPdf(Float cosThetaO) const {
    // Compute array of $A_p$ values for _cosThetaO_
    Float sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);

    // Compute $\cos \thetat$ for refracted ray
    Float sinThetaT = sinThetaO / eta;
    Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

    // Compute $\gammat$ for refracted ray
    Float etap = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
    Float sinGammaT = h / etap;
    Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));

    // Compute the transmittance _T_ of a single path through the cylinder
    Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
    std::array<Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

    // Compute $A_p$ PDF from individual $A_p$ terms
    std::array<Float, pMax + 1> apPdf;
    Float sumY =
        std::accumulate(ap.begin(), ap.end(), Float(0),
                        [](Float s, const Spectrum &ap) { return s + ap.y(); });
    for (int i = 0; i <= pMax; ++i) apPdf[i] = ap[i].y() / sumY;
    return apPdf;
}

  // (Mandy Xia) Helper function for Sample_f.
  Float HairBSDF::samplephi(int thetanum, int phionum, int phiinum, int floortheta, int floorphio, Float &pdf, Float phiunit, Float random) const{
    int num_samples =  thetanum * phionum * phiinum;
    int num_nor = phionum * phiinum;

    //std::cout<<"wavelengthindex "<<wavelengthindex<<std::endl;
    int tmpindex = wavelengthindex*num_samples + floortheta*num_nor + floorphio*phiinum;
    Float cdf[phiinum];
    for (int i = 0; i < phiinum; ++i){
      cdf[i] = CDFTABLE[tmpindex + i];
    }
    int index, indexplus;
    Float weight;
    for (index = 0; index < phiinum; ++index){
      if (random < cdf[index]){
        if (index==0)
          weight = random / cdf[0];
        else
          weight = (random-cdf[index-1])/(cdf[index]-cdf[index-1]);
        break;
      }
    }
    indexplus = Clamp(index + 1, 0, phiinum-1);
    index = Clamp(index, 0, phiinum-1);

    Float phistart = index*phiunit;
    Float phiend = indexplus*phiunit;
    Float dphi = (phiend - phistart) * weight + phistart;

    pdf = PDFTABLE[wavelengthindex*num_samples+floortheta*num_nor+floorphio*phiinum+index];

    return dphi;
  }

// (Mandy Xia) Sample using precomputation tables.
Spectrum HairBSDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u2,
                            Float *pdf, BxDFType *sampledType) const {
    // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);
    phiO += ori * M_PI/180.0;    // Take orientation into account.

    // Derive four random samples from _u2_
    Point2f u[2] = {DemuxFloat(u2[0]), DemuxFloat(u2[1])};

    // Sample $M_p$ to compute $\thetai$
    u[1][0] = std::max(u[1][0], Float(1e-5));
    Float cosTheta =
        1 + v[0] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[0]));
    Float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
    Float cosPhi = std::cos(2 * Pi * u[1][1]);
    Float sinThetaI = -cosTheta * sinThetaO + sinTheta * cosPhi * cosThetaO;
    Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));

    // Sample $N_p$ to compute $\Delta\phi$
    // Compute $\gammat$ for refracted ray
    // resolution
    int thetanum = 100;
    int phiinum = 360;
    int phionum; // incident phi

    // interpolate
    Float asino = SafeASin(std::abs(sinThetaO));
    Float unit = Pi / 2.0 / (Float)thetanum;
    int floortheta = Clamp(std::floor(asino / unit), 0, thetanum-1);
    Float phiI, dphi, pdfazimuthal, floorphio;
    Float phiOnew = phiO;
    while (phiOnew > 2*Pi) phiOnew -= 2 * Pi;
    while (phiOnew < 0) phiOnew += 2 * Pi;
    Float phiunit = M_PI * 2 / (Float)phiinum;
    if (crosssection==0){
      // circle
      phionum = 1;
      floorphio = 0;
      dphi = samplephi(thetanum, phionum, phiinum, floortheta, floorphio, pdfazimuthal, phiunit, u[0][1]);
      phiI = phiO + dphi;
    }else if (crosssection==1){
      // ellipse
      phionum = 90;
      Float phiOnewsave = phiOnew;
      if (phiOnew >= Pi/2 && phiOnew < Pi){
        phiOnew = Pi - phiOnew;
      } else if (phiOnew >= Pi && phiOnew < 3.0*Pi/2.0){
        phiOnew -= Pi;
      } else if (phiOnew >= 3.0*Pi/2.0){
        phiOnew -= Pi;
        phiOnew = Pi - phiOnew;
      }
      Float phiounit = Pi/2.0/(Float)phionum;
      floorphio = Clamp(std::floor(phiOnew / phiounit), 0, phionum-1);
      dphi = samplephi(thetanum, phionum, phiinum, floortheta, floorphio, pdfazimuthal, phiunit, u[0][1]);

      phiI = dphi;
      if (phiOnewsave >= Pi/2 && phiOnewsave < Pi){
        phiI = Pi - dphi;
      } else if (phiOnewsave >= Pi && phiOnewsave < 3.0*Pi/2.0){
        phiI += Pi;
      } else if (phiOnewsave >= 3.0*Pi/2.0){
        phiI = Pi - dphi;
        phiI += Pi;
      }
    }else{
      // nonelliptical
      phionum = 360;
      Float phiounit = Pi*2.0/(Float)phionum;
      floorphio = Clamp(std::floor(phiOnew / phiounit), 0, phionum-1);
      phiI = samplephi(thetanum, phionum, phiinum, floortheta, floorphio, pdfazimuthal, phiunit, u[0][1]);
    }
    // rotate phiI back
    phiI -= ori * M_PI/180;

    // Compute _wi_ from sampled hair scattering angles
    *wi = Vector3f(sinThetaI, cosThetaI * std::cos(phiI),
                   cosThetaI * std::sin(phiI));

    // Compute PDF for sampled hair scattering direction _wi_
    cosThetaO = std::abs(cosThetaO);
    *pdf = Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[0]) * pdfazimuthal;

    Spectrum result = f(wo, *wi);

    return result;
}

  // (Mandy Xia) Evaluate pdf using precomputation tables.
Float HairBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    // Compute hair coordinate system terms related to _wo_
    Float sinThetaO = wo.x;
    Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
    Float phiO = std::atan2(wo.z, wo.y);
    phiO += ori * M_PI/180;

    // Compute hair coordinate system terms related to _wi_
    Float sinThetaI = wi.x;
    Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
    Float phiI = std::atan2(wi.z, wi.y);
    phiI += ori * M_PI/180;

    Float phi = phiI - phiO;
    Float pdf = 0;

    int thetanum = 100;
    int phiinum = 360;
    int phionum ; // incident angle num

    Float asino = SafeASin(std::abs(sinThetaO));
    Float unit = Pi / 2.0 / (Float)thetanum;
    int floortheta = Clamp(std::floor(asino / unit), 0, thetanum-1);

    int floorphio, floorphii;
    Float phiOnew = phiO;
    while (phiOnew > 2*Pi) phiOnew -= 2 * Pi;
    while (phiOnew < 0) phiOnew += 2 * Pi;
    Float phiInew = phiI;
    while (phiInew > 2*Pi) phiInew -= 2 * Pi;
    while (phiInew < 0) phiInew += 2 * Pi;
    Float phiunit = M_PI * 2 / (Float)phiinum;
    if (crosssection==0){
      // circle
      phionum = 1;
      floorphio = 0;
      while (phi > 2*Pi) phi -= 2 * Pi;
      while (phi < 0) phi += 2 * Pi;
      floorphii = Clamp(std::floor(phi / phiunit), 0, phiinum-1);
    }else if (crosssection==1){
      // ellipse
      phionum = 90;
      //if (phiOnew < Pi/2) // don't need to change
      if (phiOnew >= Pi/2 && phiOnew < Pi){
        phiOnew = Pi - phiOnew;
        phiInew = Pi - phiInew;
        while (phiInew < 0) phiInew += 2 * Pi;
      } else if (phiOnew >= Pi && phiOnew < 3.0*Pi/2.0){
        phiOnew -= Pi;
        phiInew -= Pi;
        while (phiInew < 0) phiInew += 2 * Pi;
      } else if (phiOnew >= 3.0*Pi/2.0){
        phiOnew -= Pi;
        phiInew -= Pi;
        // treat it as second quadrant
        phiOnew = Pi - phiOnew;
        phiInew = Pi - phiInew;
        while (phiInew < 0) phiInew += 2 * Pi;
      }
      Float phiounit = Pi/2.0/(Float)phionum;
      floorphio = Clamp(std::floor(phiOnew / phiounit), 0, phionum-1);
      floorphii = Clamp(std::floor(phiInew / phiunit), 0, phiinum-1);
    }else{
      // nonelliptical
      phionum = 360;
      floorphio = Clamp(std::floor(phiOnew / phiunit), 0, phionum-1);
      floorphii = Clamp(std::floor(phiInew / phiunit), 0, phiinum-1);
    }

    unsigned long long num_samples = thetanum * phionum * phiinum;
    unsigned long long num_nor = phionum * phiinum;

    pdf = PDFTABLE[wavelengthindex*num_samples + floortheta*num_nor + floorphio*phiinum + floorphii];

    cosThetaO = std::abs(cosThetaO);
    pdf *= Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[0]);

    return pdf;
}

std::string HairBSDF::ToString() const {
    return StringPrintf(
        "[ Hair h: %f gammaO: %f eta: %f beta_m: %f beta_n: %f "
        "v[0]: %f s: %f sigma_a: ", h, gammaO, eta, beta_m, beta_n,
        v[0], s) +
        sigma_a.ToString() +
        std::string("  ]");
}

Spectrum HairBSDF::SigmaAFromConcentration(Float ce, Float cp) {
    Float sigma_a[3];
    Float eumelaninSigmaA[3] = {0.419f, 0.697f, 1.37f};
    Float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
    for (int i = 0; i < 3; ++i)
        sigma_a[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
    return Spectrum::FromRGB(sigma_a);
}

Spectrum HairBSDF::SigmaAFromReflectance(const Spectrum &c, Float beta_n) {
    Spectrum sigma_a;
    for (int i = 0; i < Spectrum::nSamples; ++i)
        sigma_a[i] = Sqr(std::log(c[i]) /
                         (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
                          10.73f * Pow<3>(beta_n) + 5.574f * Pow<4>(beta_n) +
                          0.245f * Pow<5>(beta_n)));
    return sigma_a;
}

}  // namespace pbrt
