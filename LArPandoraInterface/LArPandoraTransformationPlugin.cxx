/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraTransformationPlugin.cxx
 *
 *  @brief  Implementation of the LArPandora transformation plugin class.
 *
 *  $Log: $
 */

#include "Objects/CartesianVector.h"

#include "LArPandoraTransformationPlugin.h"

#include <cmath>

namespace lar_pandora
{

// TODO: Remove detector height from these methods (they should still work)

LArPandoraTransformationPlugin::LArPandoraTransformationPlugin(const float thetaU, const float thetaV, const float H, const float sigmaUVW,
    const float wireZPitch) :
    m_thetaU(thetaU),
    m_thetaV(thetaV),
    m_H(H),
    m_sigmaUVW(sigmaUVW),
    m_wireZPitch(wireZPitch),
    m_sinUminusV(std::sin(m_thetaU - m_thetaV)),
    m_sinUplusV(std::sin(m_thetaU + m_thetaV)),
    m_sinU(std::sin(m_thetaU)),
    m_sinV(std::sin(m_thetaV)),
    m_cosU(std::cos(m_thetaU)),
    m_cosV(std::cos(m_thetaV))
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraTransformationPlugin::~LArPandoraTransformationPlugin()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::UVtoW(const float u, const float v) const
{
    return this->UVtoZ(u, v);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::VWtoU(const float v, const float w) const
{
    return (w * m_sinUplusV - v * m_sinU + m_H * m_sinU * m_sinV) / m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::WUtoV(const float w, const float u) const
{
    return (w * m_sinUplusV - u * m_sinV + m_H * m_sinU * m_sinV) / m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::UVtoY(const float u, const float v) const
{
    return (v * m_cosU - u * m_cosV + 0.5f * m_H * m_sinUminusV) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::UVtoZ(const float u, const float v) const
{
    return (v * m_sinU + u * m_sinV - m_H * m_sinU * m_sinV) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::YZtoU(const float y, const float z) const
{
    return z * m_cosU - (y - 0.5 * m_H) * m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::YZtoV(const float y, const float z) const
{
    return z * m_cosV + (y + 0.5 * m_H) * m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::PUPVtoPW(const float pu, const float pv) const
{
    return (pu * m_sinV + pv * m_sinU) / m_sinUplusV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::PVPWtoPU(const float pv, const float pw) const
{
    return (pw * m_sinUplusV - pv * m_sinU) / m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::PWPUtoPV(const float pw, const float pu) const
{
    return (pw * m_sinUplusV - pu * m_sinV) / m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::PYPZtoPU(const float py, const float pz) const
{
    return pz * m_cosU - py * m_sinU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::PYPZtoPV(const float py, const float pz) const
{
    return pz * m_cosV + py * m_sinV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::GetSigmaUVW() const
{
    return m_sigmaUVW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraTransformationPlugin::GetWireZPitch() const
{
    return m_wireZPitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraTransformationPlugin::GetMinChiSquaredYZ(const float u, const float v, const float w, const float sigmaU, const float sigmaV,
    const float sigmaW, float &y, float &z, float &chiSquared) const
{
    // Obtain expression for chi2, differentiate wrt y and z, set both results to zero and solve simultaneously. Using Mathematica output:
    const double offsetY(0.5 * m_H * m_sinU);

    y = -((-offsetY * sigmaV * m_sinU) + (sigmaV * u * m_sinU) - (sigmaV * w * m_cosU * m_sinU) + (offsetY * sigmaW * m_cosU * m_cosV * m_sinU) -
          (sigmaW * v * m_cosU * m_cosV * m_sinU) - (offsetY * sigmaW * m_cosV * m_cosV * m_sinU) + (sigmaW * u * m_cosV * m_cosV * m_sinU) +
          (offsetY * sigmaU * m_sinV) - (sigmaU * v * m_sinV) + (offsetY * sigmaW * m_cosU * m_cosU * m_sinV) - (sigmaW * v * m_cosU * m_cosU * m_sinV) +
          (sigmaU * w * m_cosV * m_sinV) - (offsetY * sigmaW * m_cosU * m_cosV * m_sinV) + (sigmaW * u * m_cosU * m_cosV * m_sinV)) /
         ((sigmaV * m_sinU * m_sinU) + (sigmaW * m_cosV * m_cosV * m_sinU * m_sinU) + (2. * sigmaW * m_cosU * m_cosV * m_sinU * m_sinV) +
          (sigmaU * m_sinV * m_sinV) + (sigmaW * m_cosU * m_cosU * m_sinV * m_sinV));

    z = -((-sigmaV * w * m_sinU * m_sinU) + (offsetY * sigmaW * m_cosV * m_sinU * m_sinU) - (sigmaW * v * m_cosV * m_sinU * m_sinU) +
          (offsetY * sigmaW * m_cosU * m_sinU * m_sinV) - (sigmaW * v * m_cosU * m_sinU * m_sinV) + (offsetY * sigmaW * m_cosV * m_sinU * m_sinV) -
          (sigmaW * u * m_cosV * m_sinU * m_sinV) - (sigmaU * w * m_sinV * m_sinV) + (offsetY * sigmaW * m_cosU * m_sinV * m_sinV) -
          (sigmaW * u * m_cosU *m_sinV * m_sinV)) /
         ((sigmaV * m_sinU * m_sinU) + (sigmaW * m_cosV * m_cosV * m_sinU * m_sinU) + (2. * sigmaW * m_cosU * m_cosV * m_sinU * m_sinV) +
          (sigmaU * m_sinV * m_sinV) + (sigmaW * m_cosU * m_cosU * m_sinV * m_sinV));

    const float deltaU(u - LArPandoraTransformationPlugin::YZtoU(y, z));
    const float deltaV(v - LArPandoraTransformationPlugin::YZtoV(y, z));
    const float deltaW(w - z);
    chiSquared = ((deltaU * deltaU) / (sigmaU * sigmaU)) + ((deltaV * deltaV) / (sigmaV * sigmaV)) + ((deltaW * deltaW) / (sigmaW * sigmaW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraTransformationPlugin::GetProjectedYZ(const PositionAndType &hitPositionAndType, const PositionAndType &fitPositionAndType1,
    const PositionAndType &fitPositionAndType2, const float sigmaHit, const float sigmaFit, float &y, float &z, float &chiSquared) const
{
    using namespace pandora;

    const HitType hitType(hitPositionAndType.second), fitType1(fitPositionAndType1.second), fitType2(fitPositionAndType2.second);

    if (((hitType == fitType1) || (hitType == fitType2) || (fitType1 == fitType2)) ||
        ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType)) ||
        ((TPC_VIEW_U != fitType1) && (TPC_VIEW_V != fitType1) && (TPC_VIEW_W != fitType1)) ||
        ((TPC_VIEW_U != fitType2) && (TPC_VIEW_V != fitType2) && (TPC_VIEW_W != fitType2)))
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const float u((TPC_VIEW_U == hitType) ? hitPositionAndType.first : (TPC_VIEW_U == fitType1) ? fitPositionAndType1.first : fitPositionAndType2.first);
    const float v((TPC_VIEW_V == hitType) ? hitPositionAndType.first : (TPC_VIEW_V == fitType1) ? fitPositionAndType1.first : fitPositionAndType2.first);
    const float w((TPC_VIEW_W == hitType) ? hitPositionAndType.first : (TPC_VIEW_W == fitType1) ? fitPositionAndType1.first : fitPositionAndType2.first);

    const float uPrime((TPC_VIEW_U == hitType) ? LArPandoraTransformationPlugin::VWtoU(v, w) : u);
    const float vPrime((TPC_VIEW_V == hitType) ? LArPandoraTransformationPlugin::WUtoV(w, u) : v);
    const CartesianVector input(0.f, LArPandoraTransformationPlugin::UVtoY(uPrime, vPrime), LArPandoraTransformationPlugin::UVtoZ(uPrime, vPrime));

    const float position((TPC_VIEW_U == hitType) ? u : (TPC_VIEW_V == hitType) ? v : w);
    const float fitPosition((TPC_VIEW_U == hitType) ? LArPandoraTransformationPlugin::VWtoU(v, w) : (TPC_VIEW_V == hitType) ? LArPandoraTransformationPlugin::WUtoV(w, u) : LArPandoraTransformationPlugin::UVtoW(u, v));
    const CartesianVector uUnit(0., -m_cosU, m_sinU), vUnit(0., m_cosV, m_sinV), wUnit(0., 0., 1.);
    const CartesianVector &unitVector((TPC_VIEW_U == hitType) ? uUnit : (TPC_VIEW_V == hitType) ? vUnit : wUnit);

    const CartesianVector output(input + unitVector * (position - fitPosition));
    y = output.GetY();
    z = output.GetZ();

    const float sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
    const float sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
    const float sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);
    const float deltaU(u - LArPandoraTransformationPlugin::YZtoU(y, z));
    const float deltaV(v - LArPandoraTransformationPlugin::YZtoV(y, z));
    const float deltaW(w - z);
    chiSquared = ((deltaU * deltaU) / (sigmaU * sigmaU)) + ((deltaV * deltaV) / (sigmaV * sigmaV)) + ((deltaW * deltaW) / (sigmaW * sigmaW));
}

} // namespace lar_pandora
