// *********************************************************************************
// File: ParticleMpm.ipp
//
// Details: Implementation of base class for the MPM particles
//
// Dependency: MPM_ITEMS contains all template parameters
//
// Author: Samila Bandara, University of Cambridge
//
// Version: 1.0
// *********************************************************************************

// Write only the ID of this node to a given output stream
// \param[in,out]  out Output stream
// \retval         out Output stream

template<typename MPM_ITEMS>
std::ostream& mpm::ParticleMpm<MPM_ITEMS>::writeParticleId(std::ostream& out) const {
    out << id_ << " ";
    return out;
}


// Write the Coordinates of this node to a given output stream
// \param[in,out]  out Output stream
// \retval         out Output stream

template<typename MPM_ITEMS>
std::ostream& mpm::ParticleMpm<MPM_ITEMS>::writeParticleCoordinates(std::ostream& out) const {
    std::ostream_iterator<double> doubleOut(out, "  ");
    std::copy(coord_.begin(), coord_.end(), doubleOut);
    if (dim == 2) out << "  0"; // add another zero for format reasons
    out << std::endl;
    return out;
}


// Read in the coordinates from a given input stream
// \param[in,out]  inp Input stream
// \retval         inp Input stream

template<typename MPM_ITEMS>
std::istream& mpm::ParticleMpm<MPM_ITEMS>::readParticleSelf(std::istream& inp) {
    for (unsigned i = 0; i < dim; i++) {
        inp >> coord_[i];
        if (inp.peek() == ',')
            inp.ignore();
    }
    if (dim == 2) {
        double zero;
        inp >> zero;
    } // eat the third component
    return inp;
}


// Set element pointer and node pointers vector.
// \param[in]  ePtr pointer to element which has particle
// \param[in]  nVec vector containing nodes of element which has particle

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::addElementAndNodes(ElementPtr ePtr, NodeVec_ nVec) {
    ePtrOfP_ = ePtr;
    nodesOfParticle_ = nVec;
    return;
}


// Calculate shape functions and global shape function gradients only for mesh which
// lies horizontally.

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::evaluateShpfunAndGradShpfunHorMesh() {
    this->localCoordinates_();
    this->sfun_();
    // this->globalDerivatives_();     //This can also be used
    this->globalDerivativesHorMesh_();

    return;
}


// Calculate shape functions and global shape function gradients only for mesh with
// mixed types of sub meshes (horizontal or sloped).

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::evaluateShpfunAndGradShpfunMixedMesh() {
    if (meshType_ == 0) {                        // horizontal sub mesh
        this->localCoordinates_();
        this->sfun_();
        this->globalDerivativesHorMesh_();         // can use globalDerivatives_() also
    }
    else if (meshType_ == 1) {                   // sloped sub mesh
        this->localCoordinatesSlope_();
        this->sfun_();
        this->globalDerivatives_();
    }
    return;
}


// Evaluates local coordinates at particle locations for rectangular and cubic
// mesh which lies horizontally only
// \param[out]  xi_   Local cooordinate at particle position

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::localCoordinates_() {

    VecDim ctrECoord = ePtrOfP_->giveCenterCoordOfElement();
    VecDim L = ePtrOfP_->giveLengthsOfElement();

    for (unsigned i = 0; i < dim; i++) {
        xi_(i) = 2. * (coord_(i) - ctrECoord(i)) / L(i);

        if (((fabs(xi_(i)) > 0.999999) && (fabs(xi_(i)) < 1.)) || (fabs(xi_(i)) > 1.))
            xi_(i) = sign_(xi_(i));
        else if ((fabs(xi_(i)) > 0.) && (fabs(xi_(i)) < 0.000001))
            xi_(i) = 0;
    }
    return;
}


// Evaluates local coordinates at particle locations for rectangular and cubic
// mesh which lies in a slope only
// \param[out]  xi_   Local cooordinate at particle position

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::localCoordinatesSlope_() {
    VecDim ctrECoord = ePtrOfP_->giveCenterCoordOfElement();
    VecDim L = ePtrOfP_->giveLengthsOfElement();

    for (unsigned i = 0; i < dim; i++) {
        if (i < dim - 1)
            xi_(i) = 2. * (coord_(i) - ctrECoord(i)) / L(i);
        else
            xi_(i) = 2. * (coord_(i) - ctrECoord(i) - tan(slopeAngle_)
                           * (coord_(0) - ctrECoord(0))) / L(i);

        if (((fabs(xi_(i)) > 0.999999) && (fabs(xi_(i)) < 1.))
                || (fabs(xi_(i)) > 1.))
            xi_(i) = sign_(xi_(i));
        else if ((fabs(xi_(i)) > 0.) && (fabs(xi_(i)) < 0.000001))
            xi_(i) = 0.;
    }
    return;
}


// Evaluates the shape function at the given coordinate xi_ for rectangular
// and cubic mesh which lies horizontally.
// \param[in]  xi_   Local cooordinate at which the functions is evaluated
// \param[out] sfunc_  Result of all shape functions at xi

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::sfun_() {
    if (dim == 2) {
        sfunc_(0) = 0.25 * fabs((1 - xi_(0)) * (1 - xi_(1)));
        sfunc_(1) = 0.25 * fabs((1 + xi_(0)) * (1 - xi_(1)));
        sfunc_(2) = 0.25 * fabs((1 + xi_(0)) * (1 + xi_(1)));
        sfunc_(3) = 0.25 * fabs((1 - xi_(0)) * (1 + xi_(1)));
    } else if (dim == 3) {
        sfunc_(0) = fabs((1 - xi_(0)) * (1 - xi_(1)) * (1 - xi_(2))) / 8.;
        sfunc_(1) = fabs((1 + xi_(0)) * (1 - xi_(1)) * (1 - xi_(2))) / 8.;
        sfunc_(2) = fabs((1 + xi_(0)) * (1 - xi_(1)) * (1 + xi_(2))) / 8.;
        sfunc_(3) = fabs((1 - xi_(0)) * (1 - xi_(1)) * (1 + xi_(2))) / 8.;
        sfunc_(4) = fabs((1 - xi_(0)) * (1 + xi_(1)) * (1 - xi_(2))) / 8.;
        sfunc_(5) = fabs((1 + xi_(0)) * (1 + xi_(1)) * (1 - xi_(2))) / 8.;
        sfunc_(6) = fabs((1 + xi_(0)) * (1 + xi_(1)) * (1 + xi_(2))) / 8.;
        sfunc_(7) = fabs((1 - xi_(0)) * (1 + xi_(1)) * (1 + xi_(2))) / 8.;
    }
    return;
}


// The derivatives with respect to the reference coordinates evaluated at the
// given xi: \f$ \frac{d (sfunc)}{d X}(\xi) \f$. This is valid for any type of
// mesh. Use globalDerivativesHorMesh_() for general mesh to reduce time.
// \param[in]  xi       Local coordinate at which the derivates are computed
// \param[out] dsfuncDX   dsfuncDX[i,j] \f$= d (sfunc)^j / d X_i \f$

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::globalDerivatives_() {

    // shape function derivatives
    MatDimNN dsfuncDxi;
    this->sfunGrad_(dsfuncDxi);

    // get inverse of jacobi matrix
    MatDimDim Ji;
    this->jacobianInverse_(Ji);

    // global derivatives
    ublas::noalias(dsfuncDX_) = ublas::prod(Ji, dsfuncDxi);

    return;
}


// Compute the gradient w.r.t xi of the shape functions at the given xi:
// \f$ \frac{d (sfunc)}{d \xi}(\xi) \f$
// \param[in]  xi         Local coordinate at which the gradient is evaluated
// \param[out] dsfuncDxi  dsfuncDxi[i,j] \f$ = d (sfunc)^j / d \xi_i \f$

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::sfunGrad_(MatDimNN& dsfuncDxi) {

    if (dim == 2) {
        dsfuncDxi(0, 0) = -0.25 * (1 - xi_(1));
        dsfuncDxi(0, 1) =  0.25 * (1 - xi_(1));
        dsfuncDxi(0, 2) =  0.25 * (1 + xi_(1));
        dsfuncDxi(0, 3) = -0.25 * (1 + xi_(1));

        dsfuncDxi(1, 0) = -0.25 * (1 - xi_(0));
        dsfuncDxi(1, 1) = -0.25 * (1 + xi_(0));
        dsfuncDxi(1, 2) =  0.25 * (1 + xi_(0));
        dsfuncDxi(1, 3) =  0.25 * (1 - xi_(0));


    }
    else if (dim == 3) {
        dsfuncDxi(0, 0) = -0.125 * (1 - xi_(1)) * (1 - xi_(2));
        dsfuncDxi(0, 1) =  0.125 * (1 - xi_(1)) * (1 - xi_(2));
        dsfuncDxi(0, 2) =  0.125 * (1 - xi_(1)) * (1 + xi_(2));
        dsfuncDxi(0, 3) = -0.125 * (1 - xi_(1)) * (1 + xi_(2));
        dsfuncDxi(0, 4) = -0.125 * (1 + xi_(1)) * (1 - xi_(2));
        dsfuncDxi(0, 5) =  0.125 * (1 + xi_(1)) * (1 - xi_(2));
        dsfuncDxi(0, 6) =  0.125 * (1 + xi_(1)) * (1 + xi_(2));
        dsfuncDxi(0, 7) = -0.125 * (1 + xi_(1)) * (1 + xi_(2));

        dsfuncDxi(1, 0) = -0.125 * (1 - xi_(0)) * (1 - xi_(2));
        dsfuncDxi(1, 1) = -0.125 * (1 + xi_(0)) * (1 - xi_(2));
        dsfuncDxi(1, 2) = -0.125 * (1 + xi_(0)) * (1 + xi_(2));
        dsfuncDxi(1, 3) = -0.125 * (1 - xi_(0)) * (1 + xi_(2));
        dsfuncDxi(1, 4) =  0.125 * (1 - xi_(0)) * (1 - xi_(2));
        dsfuncDxi(1, 5) =  0.125 * (1 + xi_(0)) * (1 - xi_(2));
        dsfuncDxi(1, 6) =  0.125 * (1 + xi_(0)) * (1 + xi_(2));
        dsfuncDxi(1, 7) =  0.125 * (1 - xi_(0)) * (1 + xi_(2));

        dsfuncDxi(2, 0) = -0.125 * (1 - xi_(0)) * (1 - xi_(1));
        dsfuncDxi(2, 1) = -0.125 * (1 + xi_(0)) * (1 - xi_(1));
        dsfuncDxi(2, 2) =  0.125 * (1 + xi_(0)) * (1 - xi_(1));
        dsfuncDxi(2, 3) =  0.125 * (1 - xi_(0)) * (1 - xi_(1));
        dsfuncDxi(2, 4) = -0.125 * (1 - xi_(0)) * (1 + xi_(1));
        dsfuncDxi(2, 5) = -0.125 * (1 + xi_(0)) * (1 + xi_(1));
        dsfuncDxi(2, 6) =  0.125 * (1 + xi_(0)) * (1 + xi_(1));
        dsfuncDxi(2, 7) =  0.125 * (1 - xi_(0)) * (1 + xi_(1));
    }
    return;
}


// Computes inverse of Jacobi matrix for horizontal and slope mesh (2D/3D)
// \param[out]  Ji   Jacobian matrix

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::jacobianInverse_(MatDimDim& Ji) {
    VecDim L = ePtrOfP_->giveLengthsOfElement();
    Ji.clear();

    for (unsigned i = 0; i < dim; i++)
        Ji(i, i) = 2. / L(i);
    if (meshType_ == 1)
        Ji(0, dim - 1) = -2. * tan(slopeAngle_) / L(dim - 1);

    return;
}


// The derivatives with respect to the reference coordinates evaluated at the
// given xi: \f$ \frac{d (sfunc)}{d X}(\xi) \f$
// \param[in]  xi_       Local coordinate at which the derivates are computed
// \param[out] dsfuncDX_   dsfuncDX[i,j] \f$= d (sfunc)^j / d X_i \f$

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::globalDerivativesHorMesh_() {
    VecDim L = ePtrOfP_->giveLengthsOfElement();

    if (dim == 2) {
        dsfuncDX_(0, 0) = -0.5 * (1 - xi_(1)) / L(0);
        dsfuncDX_(0, 1) =  0.5 * (1 - xi_(1)) / L(0);
        dsfuncDX_(0, 2) =  0.5 * (1 + xi_(1)) / L(0);
        dsfuncDX_(0, 3) = -0.5 * (1 + xi_(1)) / L(0);

        dsfuncDX_(1, 0) = -0.5 * (1 - xi_(0)) / L(1);
        dsfuncDX_(1, 1) = -0.5 * (1 + xi_(0)) / L(1);
        dsfuncDX_(1, 2) =  0.5 * (1 + xi_(0)) / L(1);
        dsfuncDX_(1, 3) =  0.5 * (1 - xi_(0)) / L(1);


        // dsfuncDXCen_ Compute the gradient of the shape function at the center of the element
        dsfuncDXCenter_(0, 0) = -0.5 / L(0);
        dsfuncDXCenter_(0, 1) =  0.5 / L(0);
        dsfuncDXCenter_(0, 2) =  0.5 / L(0);
        dsfuncDXCenter_(0, 3) = -0.5 / L(0);

        dsfuncDXCenter_(1, 0) = -0.5 / L(1);
        dsfuncDXCenter_(1, 1) = -0.5 / L(1);
        dsfuncDXCenter_(1, 2) =  0.5 / L(1);
        dsfuncDXCenter_(1, 3) =  0.5 / L(1);
    }
    else if (dim == 3) {
        dsfuncDX_(0, 0) = -0.25 * (1 - xi_(1)) * (1 - xi_(2)) / L(0);
        dsfuncDX_(0, 1) =  0.25 * (1 - xi_(1)) * (1 - xi_(2)) / L(0);
        dsfuncDX_(0, 2) =  0.25 * (1 - xi_(1)) * (1 + xi_(2)) / L(0);
        dsfuncDX_(0, 3) = -0.25 * (1 - xi_(1)) * (1 + xi_(2)) / L(0);
        dsfuncDX_(0, 4) = -0.25 * (1 + xi_(1)) * (1 - xi_(2)) / L(0);
        dsfuncDX_(0, 5) =  0.25 * (1 + xi_(1)) * (1 - xi_(2)) / L(0);
        dsfuncDX_(0, 6) =  0.25 * (1 + xi_(1)) * (1 + xi_(2)) / L(0);
        dsfuncDX_(0, 7) = -0.25 * (1 + xi_(1)) * (1 + xi_(2)) / L(0);

        dsfuncDX_(1, 0) = -0.25 * (1 - xi_(0)) * (1 - xi_(2)) / L(1);
        dsfuncDX_(1, 1) = -0.25 * (1 + xi_(0)) * (1 - xi_(2)) / L(1);
        dsfuncDX_(1, 2) = -0.25 * (1 + xi_(0)) * (1 + xi_(2)) / L(1);
        dsfuncDX_(1, 3) = -0.25 * (1 - xi_(0)) * (1 + xi_(2)) / L(1);
        dsfuncDX_(1, 4) =  0.25 * (1 - xi_(0)) * (1 - xi_(2)) / L(1);
        dsfuncDX_(1, 5) =  0.25 * (1 + xi_(0)) * (1 - xi_(2)) / L(1);
        dsfuncDX_(1, 6) =  0.25 * (1 + xi_(0)) * (1 + xi_(2)) / L(1);
        dsfuncDX_(1, 7) =  0.25 * (1 - xi_(0)) * (1 + xi_(2)) / L(1);

        dsfuncDX_(2, 0) = -0.25 * (1 - xi_(0)) * (1 - xi_(1)) / L(2);
        dsfuncDX_(2, 1) = -0.25 * (1 + xi_(0)) * (1 - xi_(1)) / L(2);
        dsfuncDX_(2, 2) =  0.25 * (1 + xi_(0)) * (1 - xi_(1)) / L(2);
        dsfuncDX_(2, 3) =  0.25 * (1 - xi_(0)) * (1 - xi_(1)) / L(2);
        dsfuncDX_(2, 4) = -0.25 * (1 - xi_(0)) * (1 + xi_(1)) / L(2);
        dsfuncDX_(2, 5) = -0.25 * (1 + xi_(0)) * (1 + xi_(1)) / L(2);
        dsfuncDX_(2, 6) =  0.25 * (1 + xi_(0)) * (1 + xi_(1)) / L(2);
        dsfuncDX_(2, 7) =  0.25 * (1 - xi_(0)) * (1 + xi_(1)) / L(2);

        // dsfuncDXCent_ Computing the gradient of the shape fucntion at the center of the element
        dsfuncDXCenter_(0, 0) = -0.25 / L(0);
        dsfuncDXCenter_(0, 1) =  0.25 / L(0);
        dsfuncDXCenter_(0, 2) =  0.25 / L(0);
        dsfuncDXCenter_(0, 3) = -0.25 / L(0);
        dsfuncDXCenter_(0, 4) = -0.25 / L(0);
        dsfuncDXCenter_(0, 5) =  0.25 / L(0);
        dsfuncDXCenter_(0, 6) =  0.25 / L(0);
        dsfuncDXCenter_(0, 7) = -0.25 / L(0);

        dsfuncDXCenter_(1, 0) = -0.25 / L(1);
        dsfuncDXCenter_(1, 1) = -0.25 / L(1);
        dsfuncDXCenter_(1, 2) = -0.25 / L(1);
        dsfuncDXCenter_(1, 3) = -0.25 / L(1);
        dsfuncDXCenter_(1, 4) =  0.25 / L(1);
        dsfuncDXCenter_(1, 5) =  0.25 / L(1);
        dsfuncDXCenter_(1, 6) =  0.25 / L(1);
        dsfuncDXCenter_(1, 7) =  0.25 / L(1);

        dsfuncDXCenter_(2, 0) = -0.25 / L(2);
        dsfuncDXCenter_(2, 1) = -0.25 / L(2);
        dsfuncDXCenter_(2, 2) =  0.25 / L(2);
        dsfuncDXCenter_(2, 3) =  0.25 / L(2);
        dsfuncDXCenter_(2, 4) = -0.25 / L(2);
        dsfuncDXCenter_(2, 5) = -0.25 / L(2);
        dsfuncDXCenter_(2, 6) =  0.25 / L(2);
        dsfuncDXCenter_(2, 7) =  0.25 / L(2);
    }
    return;
}


// When mesh contains sloped sub  meshes this will be called.
// \param[in]  meshType  mesh type (= 1 if slope else = 0)
// \param[in]  slopeAngle  slope angle w.r.t. 0 direction (in radians)

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::addMeshTypeAndSlope(const unsigned& meshType, const double& slopeAngle) {
    meshType_ = meshType;
    slopeAngle_ = slopeAngle;
    return;
}


// send sign plus or minus
template<typename MPM_ITEMS>
double mpm::ParticleMpm<MPM_ITEMS>::sign_(double var) {
    double a;
    if (var > 0.0) a = 1.;
    if (var < 0.0) a = -1.;
    return a;
}


// clear pointers of nodes and elements
template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::clearMeshDetails()
{
    ePtrOfP_ = NULL;
    nodesOfParticle_.clear();
    meshType_ = 0;
    return;
}


// Computes the Jacobi matrix for horizontal and slope mesh (2D/3D)
// \param[out]  J   Jacobian matrix

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::jacobian_(MatDimDim& J) {
    VecDim L = ePtrOfP_->giveLengthsOfElement();
    J.clear();

    for (unsigned i = 0; i < dim; i++)
        J(i, i) = L(i) / 2.;
    if (meshType_ == 1)
        J(0, dim - 1) = L(0) * tan(slopeAngle_) / 2.;

    return;
}


// The derivatives with respect to the reference coordinates evaluated at the
// given xi: \f$ \frac{d (sfunc)}{d X}(\xi) \f$. This is valid for any type of
// mesh. But this takes 'higher computing time' than globalDerivatives_().
// \param[in]  xi       Local coordinate at which the derivates are computed
// \param[out] dsfuncDX   dsfuncDX[i,j] \f$= d (sfunc)^j / d X_i \f$

template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::globalDerivativesExtended_() {
    // nodal coordinates
    MatDimNN X;
    ePtrOfP_->nodalCoordinates(X);

    // shape function derivatives
    MatDimNN dsfuncDxi;
    this->sfunGrad_(dsfuncDxi);

    // jacobian matrix
    MatDimDim J = ublas::prod(dsfuncDxi, ublas::trans(X));

    // inverse of the jacobian
    const double detJ = corlib::inverse(J);

    // global derivatives
    ublas::noalias(dsfuncDX_) = ublas::prod(J, dsfuncDxi);

    return;
}


// Insert coordinates into an insert iterator
// \tparam IT               Type of the iterator used here
// \param[in,out]  backIter The iterator of some coordinate storage
// \retval         backIter The manipulated iterator

template<typename MPM_ITEMS>
template<typename IT>
IT mpm::ParticleMpm<MPM_ITEMS>::passParticleCoordinates(IT backIter) const {
    return std::copy(coord_.begin(), coord_.end(), backIter);
}


// store the constraints
template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::storeParticleForce(const unsigned& component, const double& value) {
    assert(component < dof);
    if (component < dof)
        particleForces_.push_back(std::make_pair(component, value));
    return;
}


// give particle forces
template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::giveParticleForces(ParticleForceVec_& GlobalParticleForces) const {
    ParticleForceVec_::const_iterator cIter = particleForces_.begin();
    ParticleForceVec_::const_iterator cEnd  = particleForces_.end();
    for (; cIter != cEnd; cIter ++) {
        const unsigned index = dofIndices_[cIter->first];
        GlobalParticleForces.push_back(std::make_pair(index,
                                       cIter->second));
    }
    return;
}


// store the constraints
template<typename MPM_ITEMS>
void mpm::ParticleMpm<MPM_ITEMS>::computeStiffness(const VecDim& xi) {
    VecNN sfunc;
    this->sfun_(xi, sfunc);
    return;
}
