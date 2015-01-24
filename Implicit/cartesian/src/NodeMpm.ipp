// *********************************************************************************
// File: NodeMpm.ipp
//
// Details: Implementation of class for nodes in MPM.
//
// Author: Krishna Kumar and Samila Bandara, University of Cambridge
//
// Version: 1.0 - Samila Bandara
// *********************************************************************************

// Constructor with id number creats nodes in NodeBasic.
// Clear the protected vectors.

template<typename MPM_ITEMS>
mpm::NodeMpm<MPM_ITEMS>::NodeMpm(const unsigned& id) : NodeBasic_(id) {
    velocitySoil_.clear();
    momentumSoil_.clear();
    accelerationSoil_.clear();
    extForceSoil_.clear();
    presForceSoil_.clear();
    intForceSoil_.clear();
    genConstraints_.clear();
    fricConstraints_.clear();
    genConstraintsSlope_.clear();
    fricConstraintsSlope_.clear();
    elementSet_.clear();
    nodeSet_.clear();
    // particleSet_.clear();

    massSoil_ = 0;

    return;
}


// intialise data and pointers to particles

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::initialise() {
    velocitySoil_.clear();
    momentumSoil_.clear();
    accelerationSoil_.clear();
    extForceSoil_.clear();
    presForceSoil_.clear();
    intForceSoil_.clear();
    // particleSet_.clear();

    massSoil_ = 0;

    return;
}


// Compute velocity using momentum balance between nodes and particles.
// This involves division by nodal mass. Then apply boundary conditions for
// velocity, where normal velocity to the boundary is zero.

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::computeSoilVelocity() {
    if (massSoil_ != 0)
        velocitySoil_ = momentumSoil_ / massSoil_;
    //    if (massSoil_ == 0)
    //      std::cout << "NOTE: m = 0 in node ' " << id_ << " ' !!! " << std::endl;

    if (genConstraints_.size())
        this->applyGeneralConstraints_(velocitySoil_);

    return;
}

//------------------------------------------------------------------------------
// Compute acceleration of soil based on momentum balance equation. Then apply
// boundary conditions (normal acceleration of the boundary = 0, tangential
// acceleration is controlled by Coulomb friction if there exists a friction
// boundary. Velocity is then calculated linearly.
// \param[in] dt   time step to calculate velocity from acceleration
// \param[in] miu  friction coefficient of boundary

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::computeSoilAccelerationAndVelocity(const double dt, const double miu) {
    if (massSoil_ != 0)
        accelerationSoil_ = (intForceSoil_ + extForceSoil_) / massSoil_;
    if (fricConstraints_.size())
        this->applyFrictionConstraints_(dt, miu);
    if (genConstraints_.size())
        this->applyGeneralConstraints_(accelerationSoil_);
    velocitySoil_ += accelerationSoil_ * dt;

    return;
}

//------------------------------------------------------------------------------
template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::updateNodalVelocity(VecDof& velocity_) {

    velocitySoil_ = velocity_;
    if (genConstraints_.size())
        this->applyGeneralConstraints_(velocitySoil_);
}

//------------------------------------------------------------------------------
// Compute velocity using momentum balance between nodes and particles.
// This involves division by nodal mass. Then apply boundary conditions for
// velocity, where normal velocity to the boundary is zero.
// If Mass of Soil is zero, setting velocity to be zero

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::computeSoilVelocityMixedMesh() {
    if (massSoil_ != 0)
        velocitySoil_ = momentumSoil_ / massSoil_;
    // if (massSoil_ == 0) {
    //      std::cout<< "NOTE: m = 0 in node ' " << id_ << " ' !!! " <<std::endl;
    // }
    if (genConstraints_.size())
        this->applyGeneralConstraints_(velocitySoil_);
    if (genConstraintsSlope_.size())
        this->applyGeneralConstraintsSlope_(velocitySoil_);

    return;
}


// Compute acceleration of soil based on momentum balance equation. Then apply
// boundary conditions (normal acceleration of the boundary = 0, tangential
// acceleration is controlled by Coulomb friction if there exists a friction
// boundary. Velocity is then calculated linearly.
// \param[in] dt   time step to calculate velocity from acceleration
// \param[in] miu  friction coefficient of boundary

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::computeSoilAccelerationAndVelocityMixedMesh(const double dt, const double miu) {
    if (massSoil_ != 0)
        accelerationSoil_ = (intForceSoil_ + extForceSoil_) / massSoil_;
    if (fricConstraints_.size())
        this->applyFrictionConstraints_(dt, miu);
    if (fricConstraintsSlope_.size())
        this->applyFrictionConstraintsSlope_(dt, miu);
    if (genConstraints_.size())
        this->applyGeneralConstraints_(accelerationSoil_);
    if (genConstraintsSlope_.size())
        this->applyGeneralConstraintsSlope_(accelerationSoil_);
    velocitySoil_ += accelerationSoil_ * dt;

    return;
}

// store coordinate direction (0/1/2) normal to boundary in a vector to apply general
//  constraints. This fix velocity and acceleration to zero in the direction specified.
//  \param[in] direction  normal direction to the boundary

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::storeGeneralConstraint(const unsigned& direction) {
    genConstraints_.push_back(direction);
    return;
}


// store coordinate direction (0/1/2) normal to boundary and sign of the normal force
// w.r.t. coordinate system (-1/+1) in a vector of std::pair to apply friction
// constraints.
// \param[in] direction  normal direction to the boundary
// \param[in] signNormDirection  sign of the normal force to the boundary
template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::storeFrictionConstraint(const unsigned& normDirection, const int& signNormDirection) {
    fricConstraints_.push_back(std::make_pair(normDirection, signNormDirection));
    return;
}


// store slope angle of the boundary in a vector to apply general constraints.
// This fix velocity and acceleration to zero normal to the boundary.
// \param[in] slopeAngle  slope angle of the submesh where node is located

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::storeGeneralConstraintSlopeBn(const double& slopeAngle) {
    genConstraintsSlope_.push_back(slopeAngle);
    return;
}


// store slope angle of the boundary and sign of the normal force w.r.t.
// slope boundary (-1/+1) in a vector of std::pair to apply friction
// constraints.
// \param[in] slopeAngle  normal direction to the boundary
// \param[in] signNormDirection  sign of the normal force to the boundary

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::storeFrictionConstraintSlopeBn(const double& slopeAngle, const int& signNormDirection) {
    fricConstraintsSlope_.push_back(std::make_pair(slopeAngle, signNormDirection));
    return;
}


// Set velocity or acceleration in the given direction in generalConstraints_
// vector to zero
// \param[in] V  velocity or acceleration

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::applyGeneralConstraints_(VecDof& V) {
    for (unsigned i = 0; i < genConstraints_.size(); i ++)
        V(genConstraints_[i]) = 0.;
    return;
}


// Limit tangential acceleration of the boundary according to Coulomb boundary law.
// frictionConstraints_ vector include the normal direction to the boundary and sign
// of the normal force into boundary (i.e. = -1 for bottom side, +1 for top side).
// tangential directions should be computed based on the normal direction and dof.
// \param[in] dt  time step
// \param[in] miu boundary friction coefficient
// \image html horBnFriction.jpg

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::applyFrictionConstraints_(const double dt, const double miu) {
    ConstraintFricVec_::const_iterator cIter = fricConstraints_.begin();

    unsigned nDir = cIter->first;             // normal direction to boundary
    int signNDir = cIter->second;             // sign of nDir wrt coord. system

    double aN, aT, vT;

    if (dim == 2) {
        unsigned tDir = (dim - 1) - nDir;     // tangential direction to boundary

        aN = accelerationSoil_(nDir);
        aT = accelerationSoil_(tDir);
        vT = velocitySoil_(tDir);

        if ((aN * signNDir) > 0.0) {
            if (vT != 0.0) {                   // kinetic friction
                double vNet = dt * aT + vT;
                double vFric = dt * miu * fabs(aN);
                if (fabs (vNet) <= vFric)
                    aT = -vT / dt;
                else
                    aT -= sign_(vNet) * miu * fabs(aN);
            } else {                            // static friction
                if (fabs(aT) <= miu * fabs(aN))
                    aT = 0.0;
                else
                    aT -= sign_(aT) * miu * fabs(aN);
            }
            accelerationSoil_(tDir) = aT;
        }

    } else if (dim == 3) {
        Mat3x2_ MatDir;
        double aTheta;
        VecDof A, V;
        unsigned tDir0, tDir1;
        double PI = std::atan(1.0) * 4;

        MatDir(0, 0) = 1;
        MatDir(0, 1) = 2;  // tangential directions for nDir = 0
        MatDir(1, 0) = 0;
        MatDir(1, 1) = 2;  // tangential directions for nDir = 1
        MatDir(2, 0) = 0;
        MatDir(2, 1) = 1;  // tangential directions for nDir = 2

        tDir0 = MatDir(nDir, 0);
        tDir1 = MatDir(nDir, 1);

        A = accelerationSoil_;
        V = velocitySoil_;

        aN = A(nDir);
        aT = sqrt(A(tDir0) * A(tDir0) + A(tDir1) * A(tDir1));
        vT = sqrt(V(tDir0) * V(tDir0) + V(tDir1) * V(tDir1));

        // if (A(tDir0) == 0) {
        if (fabs(A(tDir0)) < 0.000001) {
            if (A(tDir1) < 0.)
                aTheta = PI * (3. / 2.);
            else
                aTheta = PI / 2.;
        } else {
            aTheta = atan(A(tDir1) / A(tDir0));
        }

        if (aN * signNDir > 0.0) {
            if (vT != 0.0) {                  // kinetic friction
                ublas::bounded_vector<double, 2> vNet;
                vNet(0) = V(tDir0) + A(tDir0) * dt;  // friction applies opp to vNet
                vNet(1) = V(tDir1) + A(tDir1) * dt;
                double vTNet = sqrt(vNet(0) * vNet(0) + vNet(1) * vNet(1));
                double vTNetTheta;

                if (fabs(vNet(0)) < 0.000001) {
                    if (vNet(1) < 0)
                        vTNetTheta = PI * (3. / 2.);
                    else
                        vTNetTheta = PI / 2.;
                } else {
                    vTNetTheta = atan(vNet(1) / vNet(0));
                }

                double vFric = miu * fabs(aN) * dt;

                if (vTNet <= vFric) {
                    A(tDir0) = -V(tDir0) / dt;   // To set particle velocity to zero
                    A(tDir1) = -V(tDir1) / dt;
                } else {                         // PI gives opposite direction
                    A(tDir0) += miu * fabs(aN) * cos(PI + vTNetTheta);
                    A(tDir1) += miu * fabs(aN) * sin(PI + vTNetTheta);
                }
            } else {                              // static friction
                if (aT <= miu * fabs(aN))         // since aT is positive
                    aT = 0.0;
                else
                    aT -= miu * fabs(aN);
                A(tDir0) = aT * cos(aTheta);
                A(tDir1) = aT * sin(aTheta);
            }
        }
        accelerationSoil_ = A;
    }
    return;
}


// Set velocity or acceleration in the normal direction to the slope boundary to
// zero. The normal direction is obtained using the slope angle stored in
// genConstraintsSlope_ vector.
// \param[in] V  velocity or acceleration

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::applyGeneralConstraintsSlope_(VecDof& V) {
    double theta, Vt;
    for (unsigned i = 0; i < genConstraintsSlope_.size(); i ++) {
        theta = genConstraintsSlope_[i];
        Vt = V(0) * cos(theta) + V(dim - 1) * sin(theta);
        V(0) = Vt * cos(theta);
        V(dim - 1) = Vt * sin(theta);
    }
    return;
}


// Limit tangential acceleration of the slope boundary according to coulomb
// friction law. fricConstraintsSlope_ vector includes the slope angle and the sign
// of the normal force into boundary (i.e. = -1 for bottom side, +1 for top side).
// Normal and tangential directions are calculated using 0/1/2 direction components
// and slope angle as shown in the image.
// \param[in] dt  time step
// \param[in] miu boundary friction coefficient
// \image html slopeFriction.jpg

template<typename MPM_ITEMS>
void mpm::NodeMpm<MPM_ITEMS>::applyFrictionConstraintsSlope_(const double dt, const double miu) {
    ConstraintFricSlopeVec_::const_iterator cIter = fricConstraintsSlope_.begin();

    double beta = cIter->first;    // slope angle w.r.t. 0 direction
    int signNDir = cIter->second;  // sign of nDir wrt coord. system

    double aN, aT, vT;

    if (dim == 2) {
        aN = accelerationSoil_(1) * cos(beta) - accelerationSoil_(0) * sin(beta);
        aT = accelerationSoil_(0) * cos(beta) + accelerationSoil_(1) * sin(beta);
        vT = velocitySoil_(0) * cos(beta) + velocitySoil_(1) * sin(beta);

        if ((aN * signNDir) > 0.0) {
            if (vT != 0.0) {                                          // kinetic friction
                double vNet = dt * aT + vT;
                double vFric = dt * miu * fabs(aN);
                if (fabs (vNet) <= vFric)
                    aT = -vT / dt;
                else
                    aT -= sign_(vNet) * miu * fabs(aN);
            } else {                                                    // static friction
                if (fabs(aT) <= miu * fabs(aN))
                    aT = 0.0;
                else
                    aT -= sign_(aT) * miu * fabs(aN);
            }
            accelerationSoil_(0) = aT * cos(beta) - aN * sin(beta);
            accelerationSoil_(1) = aT * sin(beta) + aN * cos(beta);
        }

    } else if (dim == 3) {
        double aTheta, aTs, vTs;
        double PI = std::atan(1.0) * 4.;

        aN = accelerationSoil_(2) * cos(beta) - accelerationSoil_(0) * sin(beta);
        aTs = accelerationSoil_(0) * cos(beta) + accelerationSoil_(2) * sin(beta);
        aT = sqrt(aTs * aTs + accelerationSoil_(1) * accelerationSoil_(1));

        vTs = velocitySoil_(0) * cos(beta) + velocitySoil_(2) * sin(beta);
        vT = sqrt(vTs * vTs + velocitySoil_(1) * velocitySoil_(1));

        // if (A(tDir0) == 0) {
        if (fabs(aTs) < 0.000001) {
            if (accelerationSoil_(1) < 0)
                aTheta = PI * (3. / 2.);
            else
                aTheta = PI / 2.;
        } else {
            aTheta = atan(accelerationSoil_(1) / aTs);
        }

        if (aN * signNDir > 0.0) {
            if (vT != 0.0) {                               // kinetic friction
                ublas::bounded_vector< double, 3> vNet;
                for (unsigned i = 0; i < dim; i ++)
                    vNet(i) = velocitySoil_(i) + accelerationSoil_(i) * dt;
                double vTsNet = vNet(0) * cos(beta) + vNet(2) * sin(beta);
                double vTNet = sqrt(vTsNet * vTsNet + vNet(1) * vNet(1));
                double vTNetTheta;

                // if (vNet(0) == 0) {
                if (fabs(vTsNet) < 0.000001) {
                    if (vNet(1) < 0)
                        vTNetTheta = PI * (3. / 2.);
                    else
                        vTNetTheta = PI / 2.;
                } else {
                    vTNetTheta = atan(vNet(1) / vTsNet);
                }

                double vFric = miu * fabs(aN) * dt;

                if (vTNet <= vFric) {
                    accelerationSoil_(1) = -velocitySoil_(1) / dt;
                    accelerationSoil_(0) = -velocitySoil_(0) / dt;
                    accelerationSoil_(2) = -velocitySoil_(2) / dt;  // To set particle velocity to zero
                } else {                               // PI gives opposite direction
                    accelerationSoil_(1) += miu * fabs(aN) * sin(PI + vTNetTheta);
                    aTs += miu * fabs(aN) * cos(PI + vTNetTheta);
                    accelerationSoil_(0) = aTs * cos(beta);
                    accelerationSoil_(2) = aTs * sin(beta);
                }

            } else {                                   // static friction
                if (aT <= miu * fabs(aN))              // since aT is positive
                    aT = 0.0;
                else
                    aT -= miu * fabs(aN);
                accelerationSoil_(1) = aT * sin(aTheta);
                aTs = aT * cos(aTheta);
                accelerationSoil_(0) = aTs * cos(beta);
                accelerationSoil_(2) = aTs * sin(beta);
            }
        }
    }
    return;
}


// send sign plus or minus
template<typename MPM_ITEMS>
double mpm::NodeMpm<MPM_ITEMS>::sign_(double var) {
    double a;
    if (var > 0.0) a = 1.;
    if (var < 0.0) a = -1.;
    return a;
}

