 /*  NonEuclideanShell.cpp
 *  RKLibrary
 *
 *  Created by Raz Kupferman on 9/9/09.
 *  Copyright 2009 The Hebrew University. All rights reserved.
 *
 */

#include "NonEuclideanShell.H"
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <map>



/* ==============================================================================  */
/* Face Face Face Face Face Face Face Face Face Face Face Face Face Face Face Face */
/* ==============================================================================  */
/* Constructor with pointers to Nodes */
Face::Face(const TinyVector<Node*,6> &a_nodes) :
m_nodes(a_nodes),
m_faces(),
m_coordinates(),
m_edge12(),
m_edge23(),
m_edge31(),
m_nodeconnector1(),
m_nodeconnector2(),
m_nodeconnector3(),
m_nodeconnector4(),
m_nodeconnector5(),
m_nodeconnector6(),
m_Amatrix(),
m_Apivot(),
m_Bmatrix(),
m_Bpivot(),
m_abar(),
m_bbar(),
m_invabar(),
m_adjust1(1.0),
m_adjust2(1.0),
// **new members** — defaulted to zero/identity as appropriate
m_gammabar(),      // default‐constructed 2×2 of 2×2 (all zeros)
m_lambdaG(0.0),
m_muG(0.0)
{
	/* Check that the first 3 Node* are not NULL */
	assert(m_nodes(0)!=NULL && m_nodes(1)!=NULL && m_nodes(2)!=NULL);

	/* Set the Face coordinates to the average of the Nodes */
	m_coordinates = (m_nodes(0)->coordinates() + m_nodes(1)->coordinates() + m_nodes(2)->coordinates());
	m_coordinates.scale(1.0/3.0);

	/* Set the delta-Nodes */
	m_edge12      = m_nodes(1)->coordinates() - m_nodes(0)->coordinates();
	m_edge31      = m_nodes(0)->coordinates() - m_nodes(2)->coordinates();
	m_edge23      = m_nodes(2)->coordinates() - m_nodes(1)->coordinates();
}

/* ============================================================================== */
/* Assign pointers to neighboring faces */
void Face::assignNeighbors(const TinyVector<Face*,3> &a_faces)
{
	m_faces = a_faces;
}

/* ============================================================================== */
/* Initialization */
void Face::initialize(double a_thickness,
					  double a_lambda,
					  double a_mu,
					  const TinyMatrix<double,2> a_abar,
					  const TinyMatrix<double,2> a_bbar,
                      const TinyMatrix< TinyMatrix<double,2>,2 >& a_gammabar,
                      double a_lambdaG,
                      double a_muG)
{
	/* Assign the input to member data */
	m_thickness = a_thickness;
	m_lambda    = a_lambda;
	m_mu        = a_mu;
	m_abar      = a_abar;
	m_bbar      = a_bbar;
    m_gammabar  = a_gammabar;
    m_lambdaG   = a_lambdaG;
    m_muG       = a_muG;

	/* Calculate the inverse of abar */
	m_invabar  = m_abar.inverse();

	/* Calculate the (reference) area */
	m_area      = 0.5 * sqrt(m_abar.det()) * fabs(m_edge12(0)*m_edge31(1) - m_edge12(1)*m_edge31(0));

	/* Connectors to neighboring nodes */
	m_nodeconnector1 = m_nodes(0)->coordinates() - m_coordinates;
	m_nodeconnector2 = m_nodes(1)->coordinates() - m_coordinates;
	m_nodeconnector3 = m_nodes(2)->coordinates() - m_coordinates;
	m_nodeconnector4(0) = 0.023;
	m_nodeconnector4(1) = 0.321;
	m_nodeconnector5(0) = 0.432;
	m_nodeconnector5(1) = 0.923;
	m_nodeconnector6(0) = 0.754;
	m_nodeconnector6(1) = 0.147;
	if (m_nodes(3) != NULL) m_nodeconnector4 = m_nodes(3)->coordinates() - m_coordinates;
	if (m_nodes(4) != NULL) m_nodeconnector5 = m_nodes(4)->coordinates() - m_coordinates;
	if (m_nodes(5) != NULL) m_nodeconnector6 = m_nodes(5)->coordinates() - m_coordinates;

	/* Construct the A matrix */
	m_Amatrix(0,0) =   m_edge23(0)*m_edge23(0);
	m_Amatrix(0,1) = 2*m_edge23(0)*m_edge23(1);
	m_Amatrix(0,2) =   m_edge23(1)*m_edge23(1);
	m_Amatrix(1,0) =   m_edge31(0)*m_edge31(0);
	m_Amatrix(1,1) = 2*m_edge31(0)*m_edge31(1);
	m_Amatrix(1,2) =   m_edge31(1)*m_edge31(1);
	m_Amatrix(2,0) =   m_edge12(0)*m_edge12(0);
	m_Amatrix(2,1) = 2*m_edge12(0)*m_edge12(1);
	m_Amatrix(2,2) =   m_edge12(1)*m_edge12(1);
	luDecompose(m_Amatrix, m_Apivot);
	double detA = m_Amatrix.det();
	if (fabs(detA) < 1e-12) {
		std::cerr << "Warning: m_Amatrix nearly singular!" << std::endl;
}


	/* Construct the B matrix */
	m_Bmatrix(0,0) = 1.0;
	m_Bmatrix(0,1) = m_nodeconnector1(0);
	m_Bmatrix(0,2) = m_nodeconnector1(1);
	m_Bmatrix(0,3) = 0.5*m_nodeconnector1(0)*m_nodeconnector1(0);
	m_Bmatrix(0,4) = m_nodeconnector1(0)*m_nodeconnector1(1);
	m_Bmatrix(0,5) = 0.5*m_nodeconnector1(1)*m_nodeconnector1(1);
	m_Bmatrix(1,0) = 1.0;
	m_Bmatrix(1,1) = m_nodeconnector2(0);
	m_Bmatrix(1,2) = m_nodeconnector2(1);
	m_Bmatrix(1,3) = 0.5*m_nodeconnector2(0)*m_nodeconnector2(0);
	m_Bmatrix(1,4) = m_nodeconnector2(0)*m_nodeconnector2(1);
	m_Bmatrix(1,5) = 0.5*m_nodeconnector2(1)*m_nodeconnector2(1);
	m_Bmatrix(2,0) = 1.0;
	m_Bmatrix(2,1) = m_nodeconnector3(0);
	m_Bmatrix(2,2) = m_nodeconnector3(1);
	m_Bmatrix(2,3) = 0.5*m_nodeconnector3(0)*m_nodeconnector3(0);
	m_Bmatrix(2,4) = m_nodeconnector3(0)*m_nodeconnector3(1);
	m_Bmatrix(2,5) = 0.5*m_nodeconnector3(1)*m_nodeconnector3(1);
	m_Bmatrix(3,0) = 1.0;
	m_Bmatrix(3,1) = m_nodeconnector4(0);
	m_Bmatrix(3,2) = m_nodeconnector4(1);
	m_Bmatrix(3,3) = 0.5*m_nodeconnector4(0)*m_nodeconnector4(0);
	m_Bmatrix(3,4) = m_nodeconnector4(0)*m_nodeconnector4(1);
	m_Bmatrix(3,5) = 0.5*m_nodeconnector4(1)*m_nodeconnector4(1);
	m_Bmatrix(4,0) = 1.0;
	m_Bmatrix(4,1) = m_nodeconnector5(0);
	m_Bmatrix(4,2) = m_nodeconnector5(1);
	m_Bmatrix(4,3) = 0.5*m_nodeconnector5(0)*m_nodeconnector5(0);
	m_Bmatrix(4,4) = m_nodeconnector5(0)*m_nodeconnector5(1);
	m_Bmatrix(4,5) = 0.5*m_nodeconnector5(1)*m_nodeconnector5(1);
	m_Bmatrix(5,0) = 1.0;
	m_Bmatrix(5,1) = m_nodeconnector6(0);
	m_Bmatrix(5,2) = m_nodeconnector6(1);
	m_Bmatrix(5,3) = 0.5*m_nodeconnector6(0)*m_nodeconnector6(0);
	m_Bmatrix(5,4) = m_nodeconnector6(0)*m_nodeconnector6(1);
	m_Bmatrix(5,5) = 0.5*m_nodeconnector6(1)*m_nodeconnector6(1);
	luDecompose(m_Bmatrix, m_Bpivot);
}

/* ============================================================================== */
/* Calculate the position of the Face */
TinyVector<double,3> Face::position() const
{
	TinyVector<double,3> ret = (m_nodes(0)->position() + m_nodes(1)->position() + m_nodes(2)->position());
	ret.scale(1.0/3.0);
	return ret;
}

/* ============================================================================== */
/* Calculate the unit normal */
TinyVector<double,3> Face::calculateUnitNormal() const
{
	TinyVector<double,3> r1 = m_nodes(0)->position();
	TinyVector<double,3> r2 = m_nodes(1)->position();
	TinyVector<double,3> r3 = m_nodes(2)->position();
	TinyVector<double,3> ret = CrossProduct(r2-r1, r3-r2);
	ret.scale(1.0/ret.norm());
	return ret;
}
/* ============================================================================== */
/* Debugging  */



// — stretch —
double Face::stretchingEnergy() const {
    double ρ      = stretchingEnergyContentDensity();
    // double factor = std::pow(m_adjust2, -6)
    //               * m_area
    //               * m_thickness
    //               * m_adjust1;
	double invAdjust2_6 = 1.0 / (m_adjust2 * m_adjust2 * m_adjust2 *
		m_adjust2 * m_adjust2 * m_adjust2
	);

	// double base = m_thickness * m_area * m_adjust1;

	double factor = invAdjust2_6 * m_area * m_thickness * m_adjust1;
	// std::cout << "[DEBUG Face] ρ=" << ρ
    //           << " area="    << m_area
    //           << " t="       << m_thickness
    //           << " a1="      << m_adjust1
    //           << " a2="      << m_adjust2
    //           << " factor="  << factor
    //           << "\n";
    return ρ * factor;
}

// — bend —
double Face::bendingEnergy() const {
    double β      = bendingEnergyContentDensity();
    // double factor = std::pow(m_adjust2, -6)
    //               * m_area
    //               * std::pow(m_thickness * m_adjust1, 3);
	double invAdjust2_6 = 1.0 / (m_adjust2 * m_adjust2 * m_adjust2 *
		m_adjust2 * m_adjust2 * m_adjust2
	);

	double base = m_thickness  * m_adjust1;

	double factor = invAdjust2_6 * m_area * base * base * base;
    // std::cout << "[DEBUG Face] β=" << β
    //           << " area="       << m_area
    //           << " t³a1³="      << std::pow(m_thickness * m_adjust1,3)
    //           << " a2="         << m_adjust2
    //           << " factor="     << factor
    //           << "\n";
    return β * factor;
}

// — connect —
double Face::connectionEnergy() const {
    double γ      = connectionEnergyContentDensity();
    // double factor = std::pow(m_adjust2, -6)
    //               * m_area
    //               * std::pow(m_thickness * m_adjust1, 3);
	double invAdjust2_6 = 1.0 / (m_adjust2 * m_adjust2 * m_adjust2 *
		m_adjust2 * m_adjust2 * m_adjust2
	);

	double base = m_thickness  * m_adjust1;

	double factor = invAdjust2_6 * m_area * base * base * base;
    // std::cout << "[DEBUG Face] γ=" << γ
    //           << " area="        << m_area
    //           << " t³a1³="       << std::pow(m_thickness * m_adjust1,3)
    //           << " a2="          << m_adjust2
    //           << " factor="      << factor
    //           << "\n";
    return γ * factor;
}



/* ============================================================================== */
/* Calculate stretching energy density */
double Face::stretchingEnergyContentDensity() const
{
	/* Calculate the 2D metric a */
	TinyVector<double,3> r1 = m_nodes(0)->position();
	TinyVector<double,3> r2 = m_nodes(1)->position();
	TinyVector<double,3> r3 = m_nodes(2)->position();
	TinyVector<double,3> dr12 = r2 - r1;
	TinyVector<double,3> dr23 = r3 - r2;
	TinyVector<double,3> dr31 = r1 - r3;
	TinyVector<double,3> l2;
	l2(0) = innerProduct(dr23,dr23);
	l2(1) = innerProduct(dr31,dr31);
	l2(2) = innerProduct(dr12,dr12);
	TinyVector<double,3> acomp = luSolve(m_Amatrix,m_Apivot,l2);
	TinyMatrix<double,2> a;
	a(0,0) = acomp(0);
	a(0,1) = acomp(1);
	a(1,0) = acomp(1);
	a(1,1) = acomp(2);

	/* Adjust abar with adjustment parameter */
	TinyMatrix<double,2> abarAdj = m_abar;
	abarAdj.scale(m_adjust2);
	abarAdj(0,0) += 1.0 - m_adjust2;
	abarAdj(1,1) += 1.0 - m_adjust2;
	TinyMatrix<double,2> invabarAdj = abarAdj.inverse();

	/* Calculate inv(abar)(a - abar) */
	TinyMatrix<double,2> tmp  = invabarAdj*(a - abarAdj);

	// return (m_lambda*tmp.trace()*tmp.trace() + m_mu*(tmp*tmp).trace());
	double E = (m_lambda * tmp.trace() * tmp.trace() + m_mu * (tmp*tmp).trace());
// std::cout << "m_lambda = " << m_lambda
//           << ", m_mu = " << m_mu
//           << ", thickness = " << m_thickness
//           << ", area = " << m_area
//           << ", adjust1 = " << m_adjust1
//           << ", adjust2 = " << m_adjust2
//           << ", E = " << E
//           << std::endl;
	// std::cout
    //   << "[DEBUG] m_lambda=" << m_lambda
    //   << ", m_mu="     << m_mu
    //   << ", E="        << E
    //   << "\n";
	return E;
}

/* ============================================================================== */
/* Calculate bending energy density */
double Face::bendingEnergyContentDensity() const
{
	/* Save the position of the Face */
	TinyVector<double,3> my_position = position();
	TinyVector<double,3> unitnormal = calculateUnitNormal();
	TinyVector<double,6> rhs;
	for (int nodeindex=0; nodeindex<6; nodeindex++)
	{
		if (m_nodes(nodeindex) != NULL)
		{
			TinyVector<double,3> dr = m_nodes(nodeindex)->position() - my_position;
			rhs(nodeindex) = innerProduct(dr,unitnormal);
		}
		else
		{
			const TinyVector<double,2>  *ptr = NULL;
			if (nodeindex==3)      ptr = &m_nodeconnector4;
			else if (nodeindex==4) ptr = &m_nodeconnector5;
			else                   ptr = &m_nodeconnector6;

			rhs(nodeindex) = m_bbar(0,0) * (*ptr)(0) * (*ptr)(0) +
				2*m_bbar(0,1)* (*ptr)(0) * (*ptr)(1) +
				m_bbar(1,1)* (*ptr)(1) * (*ptr)(1);
		}
	}

	TinyVector<double,6> bcomp = luSolve(m_Bmatrix, m_Bpivot, rhs);
	TinyMatrix<double,2> b;
	b(0,0) = bcomp(3);
	b(0,1) = bcomp(4);
	b(1,0) = bcomp(4);
	b(1,1) = bcomp(5);

	/* Adjust abar with adjustment parameter */
	TinyMatrix<double,2> abarAdj = m_abar;
	abarAdj.scale(m_adjust2);
	abarAdj(0,0) += 1.0 - m_adjust2;
	abarAdj(1,1) += 1.0 - m_adjust2;
	TinyMatrix<double,2> invabarAdj = abarAdj.inverse();

	/* Calculate inv(abar)(b - bbar) */
	TinyMatrix<double,2> tmp  = invabarAdj*(b - m_bbar);

	return (m_lambda*tmp.trace()*tmp.trace() + m_mu*(tmp*tmp).trace()) / 3;
}
/* ============================================================================== */
/* Calculate connection energy density */
#include <iostream>    // for std::cout

double Face::connectionEnergyContentDensity() const
{
    if (m_lambdaG == 0.0 && m_muG == 0.0)
        return 0.0;

    // 1) true Christoffel
    auto Gamma = computeConnection();

    // 2) ΔΓ = Γ – Γ̄
    TinyMatrix<TinyMatrix<double,2>,2> Delta;
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          for (int l = 0; l < 2; ++l) {
            Delta(k,i)(j,l) = 
              Gamma(k,i)(j,l) - m_gammabar(k,i)(j,l);
          }
        }
      }
    }

    // 3) trace and Frobenius‐norm²
    double tr = 0.0, norm2 = 0.0;
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 2; ++i) {
        tr += Delta(k,i)(i,i);
        for (int j = 0; j < 2; ++j) {
          for (int l = 0; l < 2; ++l) {
            norm2 += Delta(k,i)(j,l) * Delta(k,i)(j,l);
          }
        }
      }
    }

    // // —— DEBUG OUTPUT ——
    // std::cout << "Face @(" << m_coordinates(0) << "," << m_coordinates(1) << "):\n"
    //           << "  trace(ΔΓ) = " << tr  << "\n"
    //           << "  ||ΔΓ||²   = " << norm2 << "\n";

    // 4) isotropic Lamé energy
    return 0.5 * m_lambdaG * (tr * tr) +   m_muG     *  norm2;
}


/* ============================================================================== */
/* Calculate the first fundamental form */
TinyVector<double,3> Face::EFG() const
{
	/* Calculate the 2D metric a */
	TinyVector<double,3> r1 = m_nodes(0)->position();
	TinyVector<double,3> r2 = m_nodes(1)->position();
	TinyVector<double,3> r3 = m_nodes(2)->position();
	TinyVector<double,3> dr12 = r2 - r1;
	TinyVector<double,3> dr23 = r3 - r2;
	TinyVector<double,3> dr31 = r1 - r3;
	TinyVector<double,3> l2;
	l2(0) = innerProduct(dr23,dr23);
	l2(1) = innerProduct(dr31,dr31);
	l2(2) = innerProduct(dr12,dr12);
	TinyVector<double,3> acomp = luSolve(m_Amatrix,m_Apivot,l2);
	// std::cout << "Face " << this << " a = " << acomp << std::endl;

	return acomp;

}


/* ============================================================================== */
/* Helper: build the 2×2 metric a from EFG */
TinyMatrix<double,2> Face::computeMetric() const
{
    // EFG() returns {E, F, G} for this face
    TinyVector<double,3> ef = EFG();
    TinyMatrix<double,2> a;
    a(0,0) = ef(0);   // E
    a(0,1) = ef(1);   // F
    a(1,0) = ef(1);   // F
    a(1,1) = ef(2);   // G
    return a;
}
/* ============================================================================== */
/* computeMetricDerivatives */
// std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
// Face::computeMetricDerivatives() const
// {
//     // ——— DEBUG: print connection stiffness params instead of da/du, da/dv ———
//     // std::cout
//     //   << "DEBUG Face @ (" << m_coordinates(0) << ", " << m_coordinates(1) << "): "
//     //   << "lambdaG = " << m_lambdaG
//     //   << ", muG = "     << m_muG
//     //   << std::endl;
//     // —————————————————————————————————————————————————————————————————————————

//     // 1) make sure neighbors exist
	
//     // assert(m_faces(0) && m_faces(1) && m_faces(2));
// 	{
//     // if any neighbor is missing, just return zeros
//     if (!m_faces(0) || !m_faces(1) || !m_faces(2)) {
// 		std::cerr << "Warning: Missing neighbor face(s) detected!" << std::endl;

//         return { TinyMatrix<double,2>(), TinyMatrix<double,2>() };
//     }
// 	}

//     // 2) pull the metric a on the three neighbor faces
//     TinyMatrix<double,2> N0 = m_faces(0)->computeMetric();
//     TinyMatrix<double,2> N1 = m_faces(1)->computeMetric();
//     TinyMatrix<double,2> N2 = m_faces(2)->computeMetric();

//     // 3) compute parametric differences Δu and Δv
//     double du = m_faces(0)->coordinates()(0)
//               - m_faces(1)->coordinates()(0);
//     double dv = m_faces(0)->coordinates()(1)
//               - m_faces(2)->coordinates()(1);
// 	if (fabs(du) < 1e-12 || fabs(dv) < 1e-12) {
// 		// handle boundary or degenerate case…
// 		return { TinyMatrix<double,2>(), TinyMatrix<double,2>() };
// 	}

//     // 4) finite‐difference the metric
//     TinyMatrix<double,2> da_du = (N0 - N1) * (1.0/du);
//     TinyMatrix<double,2> da_dv = (N0 - N2) * (1.0/dv);
// 	// std::cout << "du = " << du << ", dv = " << dv << std::endl;
// 	// if(fabs(du) < 1e-8 || fabs(dv) < 1e-8)
// 	// 	std::cerr << "Warning: Very small du/dv detected!" << std::endl;

//     return std::make_pair(da_du, da_dv);
// }

// /* computeMetricDerivatives - works well but cooses 2 neighbors*/
// std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
// Face::computeMetricDerivatives() const
// {
//     const double tol = 1e-8;

//     // 1) pull neighbor metrics (if they exist)
//     TinyMatrix<double,2> N0, N1, N2;
//     if (m_faces(0)) N0 = m_faces(0)->computeMetric();
//     if (m_faces(1)) N1 = m_faces(1)->computeMetric();
//     if (m_faces(2)) N2 = m_faces(2)->computeMetric();

//     // 2) this face’s own metric
//     TinyMatrix<double,2> a0 = computeMetric();

//     // 3) uv coords of this face and its neighbors
//     double u0  = m_coordinates(0),      v0  = m_coordinates(1);
//     double u0p = m_faces(0) ? m_faces(0)->coordinates()(0) : 0.0;
//     double u1p = m_faces(1) ? m_faces(1)->coordinates()(0) : 0.0;
//     double v0p = m_faces(0) ? m_faces(0)->coordinates()(1) : 0.0;
//     double v2p = m_faces(2) ? m_faces(2)->coordinates()(1) : 0.0;

//     // 4) allocate outputs
//     TinyMatrix<double,2> da_du, da_dv;

//     // —— ∂a/∂u using neighbors 0 & 1 —— 
//     if (m_faces(0) && m_faces(1)) {
//         double du = u0p - u1p;
//         if (fabs(du) > tol)
//             da_du = (N0 - N1) * (1.0/du);
//         else
//             da_du = TinyMatrix<double,2>();
//     }
//     else if (m_faces(0)) {
//         double du = u0p - u0;
//         if (fabs(du) > tol)
//             da_du = (N0 - a0) * (1.0/du);
//         else
//             da_du = TinyMatrix<double,2>();
//     }
//     else if (m_faces(1)) {
//         double du = u0 - u1p;
//         if (fabs(du) > tol)
//             da_du = (a0 - N1) * (1.0/du);
//         else
//             da_du = TinyMatrix<double,2>();
//     }
//     else {
//         da_du = TinyMatrix<double,2>();  // no neighbors ⇒ zero
//     }

//     // —— ∂a/∂v using neighbors 0 & 2 —— 
//     if (m_faces(0) && m_faces(2)) {
//         double dv = v0p - v2p;
//         if (fabs(dv) > tol)
//             da_dv = (N0 - N2) * (1.0/dv);
//         else
//             da_dv = TinyMatrix<double,2>();
//     }
//     else if (m_faces(0)) {
//         double dv = v0p - v0;
//         if (fabs(dv) > tol)
//             da_dv = (N0 - a0) * (1.0/dv);
//         else
//             da_dv = TinyMatrix<double,2>();
//     }
//     else if (m_faces(2)) {
//         double dv = v0 - v2p;
//         if (fabs(dv) > tol)
//             da_dv = (a0 - N2) * (1.0/dv);
//         else
//             da_dv = TinyMatrix<double,2>();
//     }
//     else {
//         da_dv = TinyMatrix<double,2>();  // no neighbors ⇒ zero
//     }

//     // 5) return the two partial derivatives
//     return { da_du, da_dv };
// }
std::pair<TinyMatrix<double,2>, TinyMatrix<double,2>>
Face::computeMetricDerivatives() const
{
    // Collect available neighbors and their UVs
    struct Neighbor {
        TinyMatrix<double,2> metric;
        double du, dv;
        double weight;
    };
    std::vector<Neighbor> neighbors;

    const TinyMatrix<double,2> a0 = computeMetric();
    const double u0 = m_coordinates(0), v0 = m_coordinates(1);

    for (int ni = 0; ni < 3; ++ni) {
        if (!m_faces(ni)) continue;
        TinyMatrix<double,2> aj = m_faces(ni)->computeMetric();
        double uj = m_faces(ni)->coordinates()(0);
        double vj = m_faces(ni)->coordinates()(1);

        double du = uj - u0;
        double dv = vj - v0;
        double dist2 = du*du + dv*dv;
        if (dist2 < 1e-12) continue;

        // Weight: inverse squared distance (safe, regularizes for closer neighbors)
        double w = 1.0 / dist2;

        neighbors.push_back({aj, du, dv, w});
    }

    // If no neighbors, return zero
    if (neighbors.empty())
        return {TinyMatrix<double,2>(), TinyMatrix<double,2>()};

    // Least-squares: Solve W*Δa = W*[du, dv] * [da_du; da_dv]
    // For each metric component, solve independently (since TinyMatrix is small)

    // Each metric is symmetric 2x2, we can unroll as [a00, a01, a11]
    // We'll fit each component separately

    TinyMatrix<double,2> da_du, da_dv;
    for (int m=0; m<2; ++m) for (int n=0; n<=m; ++n) { // (0,0), (0,1), (1,1)
        // Set up least-squares
        double Sww = 0, Swu = 0, Swv = 0, Swuu = 0, Swvv = 0, Swuv = 0;
        double Swa = 0, Swau = 0, Swav = 0;

        for (const auto& N : neighbors) {
            double da = N.metric(m,n) - a0(m,n);
            double w = N.weight;
            Sww  += w;
            Swu  += w * N.du;
            Swv  += w * N.dv;
            Swuu += w * N.du * N.du;
            Swvv += w * N.dv * N.dv;
            Swuv += w * N.du * N.dv;
            Swa  += w * da;
            Swau += w * da * N.du;
            Swav += w * da * N.dv;
        }

        // Build 2x2 system:
        // [Swuu  Swuv][X] = [Swau]
        // [Swuv  Swvv][Y]   [Swav]
        double det = Swuu*Swvv - Swuv*Swuv;
        double x = 0, y = 0;
        if (fabs(det) > 1e-14) {
            x = ( Swau*Swvv - Swav*Swuv) / det;
            y = (Swuu*Swav - Swuv*Swau) / det;
        }
        // Place into the output derivatives
        da_du(m,n) = x;
        da_dv(m,n) = y;
        if (m != n) {
            da_du(n,m) = x;
            da_dv(n,m) = y;
        }
    }
    return {da_du, da_dv};
}


/* ============================================================================== */
/* Calculate the second fundamental form */
TinyVector<double,3> Face::LMN() const
{
	/* Save the position of the Face */
	TinyVector<double,3> my_position = position();
	TinyVector<double,3> unitnormal = calculateUnitNormal();
	TinyVector<double,6> rhs;
	for (int nodeindex=0; nodeindex<6; nodeindex++)
	{
		if (m_nodes(nodeindex) != NULL)
		{
			TinyVector<double,3> dr = m_nodes(nodeindex)->position() - my_position;
			rhs(nodeindex) = innerProduct(dr,unitnormal);
		}
		else
		{
			const TinyVector<double,2>  *ptr = NULL;
			if (nodeindex==3)      ptr = &m_nodeconnector4;
			else if (nodeindex==4) ptr = &m_nodeconnector5;
			else                   ptr = &m_nodeconnector6;

			rhs(nodeindex) = m_bbar(0,0) * (*ptr)(0) * (*ptr)(0) +
				2*m_bbar(0,1)* (*ptr)(0) * (*ptr)(1) +
				m_bbar(1,1)* (*ptr)(1) * (*ptr)(1);
		}
	}

	TinyVector<double,6> bcomp = luSolve(m_Bmatrix, m_Bpivot, rhs);
	TinyVector<double,3> b;
	b(0) = bcomp(3);
	b(1) = bcomp(4);
	b(2) = bcomp(5);
	return b;
}

/* ============================================================================== */


/* ==============================================================================  */
/* Calculate the Gamma terms */
/*
   Compute the Levi–Civita connection Γᵏ_{ij} on this face from the current metric a and its derivatives.
   Uses:
     - a = computeMetric(u,v)
     - (da_du, da_dv) = computeMetricDerivatives(u,v)
   Formula:
     Γᵏ_{ij} = ½ a^{kℓ} ( ∂ᵢ a_{ℓj} + ∂ⱼ a_{ℓi} - ∂_ℓ a_{ij} )
*/
TinyMatrix<TinyMatrix<double,2>,2> Face::computeConnection() const {
	
    // 1) compute metric and its inverse
    TinyVector<double,3> ef = EFG();
    TinyMatrix<double,2> a;
    a(0,0) = ef(0);  a(0,1) = ef(1);
    a(1,0) = ef(1);  a(1,1) = ef(2);
    TinyMatrix<double,2> inva = a.inverse();

    // 2) get the two derivative‐matrices
    TinyMatrix<double,2> da_du, da_dv;
    std::tie(da_du, da_dv) = computeMetricDerivatives();

    // 3) extract inva entries into locals
    const double i00 = inva(0,0),
                 i01 = inva(0,1),
                 i11 = inva(1,1);

    // 4) alias da_du / da_dv to avoid the (i==0? …) branches inside the loop
    const TinyMatrix<double,2>* d[2] = { &da_du, &da_dv };

    TinyMatrix<TinyMatrix<double,2>,2> Gamma;

    // now only three loops: k,i,j
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 2; ++i) {
        const TinyMatrix<double,2>& di = *d[i];
        for (int j = 0; j < 2; ++j) {
          const TinyMatrix<double,2>& dj = *d[j];

          // unroll l=0 and l=1 by hand:
          //   term(l) = ∂ᵢa_{ℓj} + ∂ⱼa_{ℓi} – ∂_ℓa_{ij}
          double term0 = di(0,j) + dj(0,i) - (*d[0])(i,j);
          double term1 = di(1,j) + dj(1,i) - (*d[1])(i,j);

          // sum = ∑ₗ a^{kℓ} * term(l)
          double sum = (k==0 ? (i00*term0 + i01*term1)
                             : (i01*term0 + i11*term1));

          Gamma(k,i)(0,j) = 0.5 * sum;
        }
      }
    }

    return Gamma;
}
// Note: the original code had a commented‐out version of this function
// that used a different indexing scheme for Gamma, which was slow.
// The current version uses the correct indexing as per the formula above.
// TinyMatrix<TinyMatrix<double,2>,2> Face::computeConnection() const {
//     TinyVector<double,3> ef = EFG();  // First fundamental form

//     TinyMatrix<double,2> a;
//     a(0,0) = ef(0);  a(0,1) = ef(1);
//     a(1,0) = ef(1);  a(1,1) = ef(2);
//     TinyMatrix<double,2> inva = a.inverse();

//     TinyMatrix<double,2> da_du, da_dv;
//     std::tie(da_du, da_dv) = computeMetricDerivatives();  // Per face

//     TinyMatrix<TinyMatrix<double,2>,2> Gamma;

//     for (int k = 0; k < 2; ++k) {
//         for (int i = 0; i < 2; ++i) {
//             for (int j = 0; j < 2; ++j) {
//                 double sum = 0.0;
//                 for (int l = 0; l < 2; ++l) {
//                     const TinyMatrix<double,2>& d_i = (i == 0 ? da_du : da_dv);
//                     const TinyMatrix<double,2>& d_j = (j == 0 ? da_du : da_dv);
//                     const TinyMatrix<double,2>& d_l = (l == 0 ? da_du : da_dv);

//                     double term = d_i(l,j) + d_j(l,i) - d_l(i,j);
//                     sum += inva(k,l) * term;
//                 }
//                 Gamma(k,i)(0,j) = 0.5 * sum;  // ✅ specify both indices
// 				// Gamma(k,i)(j,0) = 0.5 * sum; 

//             }
//         }
//     }

//     return Gamma;
// }
/* ============================================================================== */
/* Calculate the energy gradient */
void Face::setForce() {
    double ep = 1.e-6;
    for (int i=0; i<6; i++) {
        if (m_nodes(i) != NULL) {
            for (int comp=0; comp<3; comp++) {
                m_nodes(i)->position(comp) += ep;
                double Eplus = energy();
                m_nodes(i)->position(comp) -= 2*ep;
                double Eminus = energy();
                double grad = 0.5*(Eplus-Eminus)/ep;

                // Diagnostic print statement here:
                // std::cout << "Gradient component [" << i << "][" << comp << "] = " << grad << std::endl;

                m_nodes(i)->force(comp) += grad;
                m_nodes(i)->position(comp) += ep;
            }
        }
    }
}

/* ============================================================================== */
/* NonEuclideanShell NonEuclideanShell NonEuclideanShell NonEuclideanShell    */
/* ============================================================================== */
/* Constructor */
NonEuclideanShell::NonEuclideanShell(const std::string &a_nodesFileName,
									 const std::string &a_facesFileName) :
m_nodes(),
m_faces(),
m_verbosity(2)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::NonEuclideanShell()");

	int		numberNodes, numberFaces;
	int		n1, n2, n3, n4, n5, n6, f1, f2, f3;
	int		index;
	double	u,v;
	int		isfixed;

	/* Take care of Nodes */
	TextFileHandle nodesFileHandle(a_nodesFileName,FileHandle::OPEN_RD);
	nodesFileHandle.read(numberNodes);
	m_nodes = Vector<Node*>(numberNodes);

	if (m_verbosity>1)
		std::cout << "NonEuclideanShell::NonEuclideanShell()   Number of Nodes = " << numberNodes << std::endl;
    // ——— DEBUG: collect counts of each isfixed value ———
    std::map<int,int> flagCounts;

	for (int i=0; i<numberNodes; i++)
	{
		nodesFileHandle.read(index);
		nodesFileHandle.read(u);
		nodesFileHandle.read(v);
		nodesFileHandle.read(isfixed);
		m_nodes(i)  = new Node();
        // ——— DEBUG: print each one as it arrives ———
        // std::cerr << "[DEBUG] Read vertex " << index
        //           << " → u=" << u << ", v=" << v
        //           << ", isfixed=" << isfixed << "\n";
        // tally for summary
        flagCounts[isfixed]++;
		m_nodes(i)->coordinates(0) = u;
		m_nodes(i)->coordinates(1) = v;
		m_nodes(i)->fixed() = isfixed;
	}
	nodesFileHandle.close();
	  // ——— DEBUG: print a summary of all flags seen ———
    std::cerr << "[DEBUG] isfixed value summary:\n";
    for (auto &kv : flagCounts) {
        std::cerr << "   " << kv.first << " → " << kv.second << " vertices\n";
    }

	/* Take care of Faces */
	TextFileHandle facesFileHandle(a_facesFileName,FileHandle::OPEN_RD);
	facesFileHandle.read(numberFaces);
	m_faces = Vector<Face*>(numberFaces);

	if (m_verbosity>1)
		std::cout << "NonEuclideanShell::NonEuclideanShell()   Number of Faces = " << numberFaces << std::endl;

	for (int i=0; i<numberFaces; i++)
	{
		facesFileHandle.read(index);
		facesFileHandle.read(n1);
		facesFileHandle.read(n2);
		facesFileHandle.read(n3);
		facesFileHandle.read(n4);
		facesFileHandle.read(n5);
		facesFileHandle.read(n6);
		facesFileHandle.read(f1);
		facesFileHandle.read(f2);
		facesFileHandle.read(f3);

		TinyVector<Node*, 6> nodeVec;
		nodeVec(0) = m_nodes(n1);
		nodeVec(1) = m_nodes(n2);
		nodeVec(2) = m_nodes(n3);
		nodeVec(3) = (n4==-1) ? NULL : m_nodes(n4);
		nodeVec(4) = (n5==-1) ? NULL : m_nodes(n5);
		nodeVec(5) = (n6==-1) ? NULL : m_nodes(n6);

		m_faces(i) = new Face(nodeVec);
	}
	facesFileHandle.close();

	/* Take care of Faces */
	TextFileHandle facesFileHandle2(a_facesFileName,FileHandle::OPEN_RD);
	facesFileHandle2.read(numberFaces);

	for (int i=0; i<numberFaces; i++)
	{
		facesFileHandle2.read(index);
		facesFileHandle2.read(n1);
		facesFileHandle2.read(n2);
		facesFileHandle2.read(n3);
		facesFileHandle2.read(n4);
		facesFileHandle2.read(n5);
		facesFileHandle2.read(n6);
		facesFileHandle2.read(f1);
		facesFileHandle2.read(f2);
		facesFileHandle2.read(f3);

		TinyVector<Face*, 3> faceVec;
		faceVec(0) = (f1==-1) ? NULL : m_faces(f1);
		faceVec(1) = (f2==-1) ? NULL : m_faces(f2);
		faceVec(2) = (f3==-1) ? NULL : m_faces(f3);

		m_faces(i)->assignNeighbors(faceVec);
	}
	facesFileHandle2.close();

}

/* ============================================================================== */
/* Destructor */
NonEuclideanShell::~NonEuclideanShell()
{
	for (int i=0; i<m_nodes.length(); i++) delete m_nodes(i);
	for (int i=0; i<m_faces.length(); i++) delete m_faces(i);
}

/* ============================================================================== */
/* Set reference forms and parameters */
void NonEuclideanShell::setParameters(double (*f_thickness)(double, double),
									  double (*f_lambda)(double, double),
									  double (*f_mu)(double, double),
									  TinyMatrix<double,2> (*f_abar)(double, double),
									  TinyMatrix<double,2> (*f_bbar)(double, double),
									  TinyVector<double,3> (*f_Pos0)(double, double),
									  TinyMatrix<TinyMatrix<double,2>,2> (*f_gammabar)(double, double),
                                      double lambdaG,
                                      double muG)
{
	for (int i=0; i<m_faces.length(); i++)
	{
		double u         = m_faces(i)->coordinates(0);
		double v         = m_faces(i)->coordinates(1);
		double thickness = f_thickness(u,v);
		double lambda    = f_lambda(u,v);
		double mu        = f_mu(u,v);
		TinyMatrix<double,2> abar      = f_abar(u,v);
		TinyMatrix<double,2> bbar      = f_bbar(u,v);
        TinyMatrix<TinyMatrix<double,2>,2> gammabar = f_gammabar(u, v);

		m_faces(i)->initialize(thickness, lambda, mu, abar, bbar, gammabar, lambdaG, muG);
    }

    // --- initialize node positions (unchanged) ---
    for (int i = 0; i < m_nodes.length(); ++i) 
	{
        double u = m_nodes(i)->coordinates()(0);
        double v = m_nodes(i)->coordinates()(1);
        // auto   P = f_Pos0(u, v);
        // m_nodes(i)->position()(0) = P(0);
        // m_nodes(i)->position()(1) = P(1);
        // m_nodes(i)->position()(2) = P(2);
		m_nodes(i)->position(0) = f_Pos0(u,v)(0);
		m_nodes(i)->position(1) = f_Pos0(u,v)(1);
		m_nodes(i)->position(2) = f_Pos0(u,v)(2);
    }

    // --- fix offsets (unchanged) ---
    for (int i = 0; i < m_nodes.length(); ++i) {
        int j = m_nodes(i)->fixed();
        if (j >= 0) {
            m_nodes(i)->offset() =
                m_nodes(i)->position() - m_nodes(j)->position();
        }
    }
}
/* ============================================================================== */
/* set the adjustment parameters */
void NonEuclideanShell::setAdjust(double a_adjust1, double a_adjust2)
{
	for (int i=0; i<m_faces.length(); i++)
	{
		m_faces(i)->setAdjust(a_adjust1, a_adjust2);
	}
}

/* ============================================================================== */
/* Initialization: set initial configuration */
void NonEuclideanShell::defaultInitialization()
{
	for (int i=0; i<m_nodes.length(); i++)
	{
		m_nodes(i)->position(0) = m_nodes(i)->coordinates(0);
		m_nodes(i)->position(1) = m_nodes(i)->coordinates(1);
		m_nodes(i)->position(2) = 0.1*sin(m_nodes(i)->coordinates(0));
	}
}

/* ============================================================================== */
/* Restart from file */
void NonEuclideanShell::restart(BinaryFileHandle* a_fh)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::restart()");

	for (int i=0; i<m_nodes.length(); i++)
	{
		float X,Y,Z;
		a_fh->read(&X);
		a_fh->read(&Y);
		a_fh->read(&Z);
		m_nodes(i)->position(0) = X;
		m_nodes(i)->position(1) = Y;
		m_nodes(i)->position(2) = Z;
	}
}

/* ============================================================================== */
/* Initialization: set forces to zero */
void NonEuclideanShell::initializeForce()
{
	for (int i=0; i<m_nodes.length(); i++)
		m_nodes(i)->force().setToZero();
}

/* ============================================================================== */
/* Calculate the stretching energy */
// double NonEuclideanShell::stretchingEnergy() const
// {
// 	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::stretchingEnergy()");
	
// 	double energy=0;
// 	for (int i=0; i<m_faces.length(); i++) {
// 		energy += m_faces(i)->stretchingEnergy();
// 		// std::cout << "Face " << i << " stretch energy: " << m_faces(i)->stretchingEnergy() << std::endl;
// 	}

// 	if (m_verbosity>1)
// 		std::cout << "NonEuclideanShell::stretchingEnergy()    = " << energy << std::endl;
// 	// std::cout << "NonEuclideanShell::stretchingEnergy()    = " << energy << std::endl;

// 	return energy;
// }
double NonEuclideanShell::stretchingEnergy() const {
    double energy = 0;
    for (int i = 0; i < m_faces.length(); ++i) {
        // double density = m_faces(i)->stretchingEnergyContentDensity();
        // double A       = m_faces(i)->area();
        double Eface   = m_faces(i)->stretchingEnergy();
        // std::cout << "[DEBUG] Face " << i
        //           << " ρ=" << density
        //           << "  A=" << A
        //           << " E=" << Eface
        //           << std::endl;
        energy += Eface;
    }
    // std::cout << "[DEBUG] Total stretching = " << energy << std::endl;
    return energy;
}


/* ============================================================================== */
/* Calculate the bending energy */
double NonEuclideanShell::bendingEnergy() const
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::bendingEnergy()");

	double energy=0;
	for (int i=0; i<m_faces.length(); i++)
		energy += m_faces(i)->bendingEnergy();

	if (m_verbosity>1)
		std::cout << "NonEuclideanShell::bendingEnergy()    = " << energy << std::endl;

	return energy;
}
/* ============================================================================== */
/* Calculate the connection (γ) energy */
double NonEuclideanShell::connectionEnergy() const
{
    if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::connectionEnergy()");

    double energy = 0.0;
    for (int i = 0; i < m_faces.length(); ++i)
        energy += m_faces(i)->connectionEnergy();

    if (m_verbosity>1)
        std::cout << "NonEuclideanShell::connectionEnergy()    = " << energy << std::endl;

    return energy;
}

/* ============================================================================== */
/* Calculate the total energy */
double NonEuclideanShell::energy() const
{
    double energy;
    if (m_faces.length() > 0 && m_faces(0)->getLambdaG() == 0.0 && m_faces(0)->getMuG() == 0.0) {
    // std::cout << "[DEBUG] Skipping connection energy (lambdaG and muG are zero for all faces)\n";
    energy = stretchingEnergy() + bendingEnergy();
	} else {
    energy = stretchingEnergy() + bendingEnergy() + connectionEnergy();
	}


    if (m_verbosity > 1)
        std::cout << "NonEuclideanShell::energy()    = " << energy << std::endl;

    return energy;
}

/* ============================================================================== */
/* Calculate size of optimization problem */
int NonEuclideanShell::SizeOfOptimizationProblem() const
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::SizeOfOptimizationProblem()");

	int sop = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2) sop+=3;
	}

	return sop;
}

/* ============================================================================== */
/* Calculate the stretching forces */
void NonEuclideanShell::setForce()
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::setForce()");

	for (int i=0; i<m_faces.length(); i++)
		m_faces(i)->setForce();

}

/* ============================================================================== */
/* Dump results to file (binary) */
void NonEuclideanShell::DumpStateBinaryFormat(BinaryFileHandle* a_fh)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::DumpStateBinaryFormat()");

	for (int i=0; i<m_nodes.length(); i++)
	{
		float X = (float)m_nodes(i)->position(0);
		float Y = (float)m_nodes(i)->position(1);
		float Z = (float)m_nodes(i)->position(2);
		a_fh->write(&X);
		a_fh->write(&Y);
		a_fh->write(&Z);
	}
}
/* ============================================================================== */
/* Dump results to file (text) */
void NonEuclideanShell::DumpStateTextFormat(TextFileHandle* a_node, TextFileHandle* a_face, int counter)
{
	if (m_verbosity>3) Errors::StepIn("NonEuclideanShell::DumpStateTextFormat()");

	for (int i=0; i<m_nodes.length(); i++)
	{
		double X = m_nodes(i)->position(0);
		double Y = m_nodes(i)->position(1);
		double Z = m_nodes(i)->position(2);
		a_node->write(i);
		a_node->tab();
		a_node->write(X);
		a_node->tab();
		a_node->write(Y);
		a_node->tab();
		a_node->write(Z);
		a_node->tab();
		a_node->write(counter);
		a_node->newline();
	}

	for (int i=0; i<m_faces.length(); i++)
	{
		double Es = m_faces(i)->stretchingEnergyContentDensity();
		double Eb = m_faces(i)->bendingEnergyContentDensity();
		double Eg = m_faces(i)->connectionEnergyContentDensity();
		a_face->write(i);
		a_face->tab();
		a_face->write(Es);
		a_face->tab();
		a_face->write(Eb);
		a_face->tab();
		a_face->write(Eg);            // ← write the γ‐energy
        a_face->tab();
		a_face->write(counter);
		a_face->newline();
	}

}
/* ============================================================================== */
/* Dump forms to file (text) */
// void NonEuclideanShell::DumpFormsTextFormat(TextFileHandle* a_face)
// {
// 	for (int i=0; i<m_faces.length(); i++)
// 	{
// 		TinyVector<double,3> EFGv = m_faces(i)->EFG();
// 		TinyVector<double,3> LMNv = m_faces(i)->LMN();
// 		TinyVector<double,2> UVv = m_faces(i)->coordinates();
// 		a_face->write(i);
// 		a_face->tab();
// 		a_face->write(UVv(0));
// 		a_face->tab();
// 		a_face->write(UVv(1));
// 		a_face->tab();
// 		a_face->write(EFGv(0));
// 		a_face->tab();
// 		a_face->write(EFGv(1));
// 		a_face->tab();
// 		a_face->write(EFGv(2));
// 		a_face->tab();
// 		a_face->write(LMNv(0));
// 		a_face->tab();
// 		a_face->write(LMNv(1));
// 		a_face->tab();
// 		a_face->write(LMNv(2));
// 		a_face->newline();
// 	}
// }
/* ============================================================================== */
/* Dump forms (E, F, G, L, M, N) + Γ⁽ᵏ⁾ᵢⱼ (2×2×2 = 8 values) to file (text) */
void NonEuclideanShell::DumpFormsTextFormat(TextFileHandle* a_face)
{
    for (int i = 0; i < m_faces.length(); ++i)
    {
        auto   EFGv   = m_faces(i)->EFG();
        auto   LMNv   = m_faces(i)->LMN();
        auto   UVv    = m_faces(i)->coordinates();
        double Es  = m_faces(i)->stretchingEnergyContentDensity();
		double Eb  = m_faces(i)->bendingEnergyContentDensity();
		double Eg  = m_faces(i)->connectionEnergyContentDensity();

        a_face->write(i);        a_face->tab();
        a_face->write(UVv(0));   a_face->tab();
        a_face->write(UVv(1));   a_face->tab();

        a_face->write(EFGv(0));  a_face->tab();
        a_face->write(EFGv(1));  a_face->tab();
        a_face->write(EFGv(2));  a_face->tab();

        a_face->write(LMNv(0));  a_face->tab();
        a_face->write(LMNv(1));  a_face->tab();
        a_face->write(LMNv(2));  a_face->tab();

        a_face->write(Es);  a_face->tab();
		a_face->write(Eb);  a_face->tab();
		a_face->write(Eg);  
        a_face->newline();
    }
}
/* ============================================================================== */
/* I/O of state and force. Needed for external optimization procedure */
void NonEuclideanShell::setPositionVector(const double *ptr)
{
	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			m_nodes(i)->position(0) = ptr[3*j];
			m_nodes(i)->position(1) = ptr[3*j+1];
			m_nodes(i)->position(2) = ptr[3*j+2];
			j += 1;
		}
	}

	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(i)->position() = m_nodes(j)->position() + m_nodes(i)->offset();
		}
	}
}

void NonEuclideanShell::setPositionVector(const gsl_vector *a_vec)
{
	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			m_nodes(i)->position(0) = gsl_vector_get(a_vec,3*j);
			m_nodes(i)->position(1) = gsl_vector_get(a_vec,3*j+1);
			m_nodes(i)->position(2) = gsl_vector_get(a_vec,3*j+2);
			j += 1;
		}
	}

	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(i)->position() = m_nodes(j)->position() + m_nodes(i)->offset();
		}
	}
}
/* ============================================================================== */
/* I/O of state and force. Needed for external optimization procedure */
void NonEuclideanShell::getPositionVector(gsl_vector *a_vec)
{
	int j = 0;
	double* ptr = gsl_vector_ptr(a_vec, 0);
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			ptr[3*j]   = m_nodes(i)->position(0);
			ptr[3*j+1] = m_nodes(i)->position(1);
			ptr[3*j+2] = m_nodes(i)->position(2);
			j += 1;
		}
	}
}
/* ============================================================================== */
/* return the energy given an array containing the position */
double NonEuclideanShell::getEnergy(const gsl_vector *a_vec)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergy()");

	const double* ptr = gsl_vector_const_ptr(a_vec, 0);
	setPositionVector(ptr);
	return energy();
}
/* ============================================================================== */
/* return the energy gradient given an array containing the position */
void NonEuclideanShell::getEnergyGradient(const gsl_vector *a_state, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyGradient()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	initializeForce();
	setForce();

	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(j)->force() += m_nodes(i)->force();
		}
	}

	j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			ret_ptr[3*j]   = m_nodes(i)->force(0);
			ret_ptr[3*j+1] = m_nodes(i)->force(1);
			ret_ptr[3*j+2] = m_nodes(i)->force(2);
			j += 1;
		}
	}
}
/* ============================================================================== */
/* return the energy and its gradient given an array containing the position */
void NonEuclideanShell::getEnergyAndEnergyGradient(const gsl_vector *a_state,
								  double *a_energy, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyAndEnergyGradient()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	initializeForce();
	setForce();

	int j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		j = m_nodes(i)->fixed();
		if (j >= 0)
		{
			m_nodes(j)->force() += m_nodes(i)->force();
		}
	}

	j = 0;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			ret_ptr[3*j]   = m_nodes(i)->force(0);
			ret_ptr[3*j+1] = m_nodes(i)->force(1);
			ret_ptr[3*j+2] = m_nodes(i)->force(2);
			j += 1;
		}
	}
	*a_energy = energy();
}
/* ============================================================================== */
/* return the energy gradient given an array containing the position, Full calculation */
void NonEuclideanShell::getEnergyGradientFull(const gsl_vector *a_state, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyGradientFull()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	double ep    = 1.e-6;
	double Einit = energy();
	int j = 0;

	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			for (int comp=0; comp<3; comp++)
			{
				m_nodes(i)->position(comp) += ep;
				double Eplus = energy();
				m_nodes(i)->position(comp) -= 2*ep;
				double Eminus = energy();
				ret_ptr[3*j+comp] = 0.5*(Eplus-Eminus)/ep;
				m_nodes(i)->position(comp) += ep;
			}
			j += 1;
		}
	}
}
/* ============================================================================== */
/* return the energy and its gradient given an array containing the position, Full calculation */
void NonEuclideanShell::getEnergyAndEnergyGradientFull(const gsl_vector *a_state,
								  double *a_energy, gsl_vector *a_gradient)
{
	if (m_verbosity>3)  Errors::StepIn("NonEuclideanShell::getEnergyAndEnergyGradientFull()");

	const double* ptr     = gsl_vector_const_ptr(a_state, 0);
	double*       ret_ptr = gsl_vector_ptr(a_gradient, 0);
	setPositionVector(ptr);
	double ep    = 1.e-6;
	double Einit = energy();
	int j = 0;

	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			for (int comp=0; comp<3; comp++)
			{
				m_nodes(i)->position(comp) += ep;
				double Eplus = energy();
				m_nodes(i)->position(comp) -= 2*ep;
				double Eminus = energy();
				ret_ptr[3*j+comp] = 0.5*(Eplus-Eminus)/ep;
				m_nodes(i)->position(comp) += ep;
			}
			j += 1;
		}
	}
	*a_energy = energy();
}

/* ============================================================================== */
/* Test the gradient calculation                                                 */
void NonEuclideanShell::testGradient()
{
	int    noDOF = SizeOfOptimizationProblem();
	double grad[noDOF], gradFull[noDOF];

	/* calculate force */
	initializeForce();
	setForce();
	int j = 0;

	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			grad[3*j]   = m_nodes(i)->force(0);
			grad[3*j+1] = m_nodes(i)->force(1);
			grad[3*j+2] = m_nodes(i)->force(2);
			j += 1;
		}
	}

	j = 0;

	/* calculate force directly */
	double ep    = 1.e-6;
	for (int i=0; i<m_nodes.length(); i++)
	{
		if (m_nodes(i)->fixed()==-2)
		{
			for (int comp=0; comp<3; comp++)
			{
				m_nodes(i)->position(comp) += ep;
				double Eplus = energy();
				m_nodes(i)->position(comp) -= 2*ep;
				double Eminus = energy();
				gradFull[3*j+comp] = 0.5*(Eplus-Eminus)/ep;
				m_nodes(i)->position(comp) += ep;
			}
			j += 1;
		}
	}

	/* output */
	for (int i=0; i<noDOF; i++)
	{
		// std::cout << grad[i] << " " << gradFull[i] << " " <<  fabs(gradFull[i]-grad[i]) << std::endl;
		std::cout << grad[i] << " " << gradFull[i] << " " << fabs(gradFull[i]-grad[i])
          << " " << (fabs(gradFull[i]) > 1e-10 ? fabs(gradFull[i]-grad[i])/fabs(gradFull[i]) : 0.0)
          << std::endl;	

	}


	}
	const Vector<Face*>& NonEuclideanShell::getFaces() const {
		return m_faces;
	}
/* ============================================================================== */
void NonEuclideanShell::checkNodePositions() const {
    for (int i=0; i<m_nodes.length(); ++i) {
        for (int d=0; d<3; ++d) {
            if (!std::isfinite(m_nodes(i)->position(d))) {
                std::cerr << "Invalid node " << i << " pos: "
                          << m_nodes(i)->position(0) << ", "
                          << m_nodes(i)->position(1) << ", "
                          << m_nodes(i)->position(2) << std::endl;
                exit(1);
            }
        }
    }
}

