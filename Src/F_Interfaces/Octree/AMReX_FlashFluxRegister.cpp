#include <AMReX_FlashFluxRegister.H>

namespace amrex {

FlashFluxRegister::FlashFluxRegister (const BoxArray& fba, const BoxArray& cba,
                                      const DistributionMapping& fdm, const DistributionMapping& cdm,
                                      const Geometry& fgeom, const Geometry& cgeom,
                                      IntVect const& ref_ratio, int nvar)
{
    define(fba,cba,fdm,cdm,fgeom,cgeom,ref_ratio,nvar);
}

void FlashFluxRegister::define (const BoxArray& fba, const BoxArray& cba,
                                const DistributionMapping& fdm, const DistributionMapping& cdm,
                                const Geometry& fgeom, const Geometry& cgeom,
                                IntVect const& a_ref_ratio, int nvar)
{
    constexpr int ref_ratio = 2;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_ref_ratio == IntVect(ref_ratio),
                                     "FlashFluxRegister: refinement ratio != 2");

    m_fine_grids = fba;
    m_crse_grids = cba;
    m_fine_dmap = fdm;
    m_crse_dmap = cdm;
    m_fine_geom = fgeom;
    m_crse_geom = cgeom;
    m_ncomp = nvar;

    const int myproc = ParallelDescriptor::MyProc();

    // For a fine Box, there is at most one face per direction abutting coarse level.
    {
        BoxArray const& fndba = amrex::convert(fba,IntVect::TheNodeVector());
        Array<BoxList,AMREX_SPACEDIM> bl;
        Array<Vector<int>,AMREX_SPACEDIM> procmap;
        Array<Vector<int>,AMREX_SPACEDIM> my_global_indices;
        std::vector<std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = m_fine_geom.periodicity().shiftIntVect();
        Vector<std::pair<Orientation,Box> > faces;
        for (int i = 0, nfines = fba.size(); i < nfines; ++i)
        {
            Box const& ccbx = fba[i];
            faces.clear();
            for (OrientationIter orit; orit.isValid(); ++orit) {
                const Orientation ori = orit();
                faces.emplace_back(std::make_pair(ori,amrex::bdryNode(ccbx,ori)));
            }
            Box const& ndbx = amrex::surroundingNodes(ccbx);
            for (auto const& shift : pshifts) {
                if (!faces.empty()) {
                    fndba.intersections(ndbx+shift, isects);
                    for (auto const& isect : isects) {
                        Box const& b = isect.second-shift;
                        faces.erase(std::remove_if(faces.begin(), faces.end(),
                                                   [&] (std::pair<Orientation,Box> const& x)
                                                   { return x.second == b; }),
                                    faces.end());
                    }
                }
            }
            for (auto const& face : faces) {  // coarse/fine boundary faces
                const int dir = face.first.coordDir();
                bl[dir].push_back(amrex::coarsen(face.second,ref_ratio));
                procmap[dir].push_back(fdm[i]);
                if (fdm[i] == myproc) my_global_indices[dir].push_back(i);
            }
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (bl[idim].isNotEmpty()) {
                BoxArray faces_ba(std::move(bl[idim]));
                DistributionMapping faces_dm(std::move(procmap[idim]));
                m_fine_fluxes[idim].define(faces_ba, faces_dm, m_ncomp, 0);
                for (MFIter mfi(m_fine_fluxes[idim]); mfi.isValid(); ++mfi) {
                    const int gi = my_global_indices[idim][mfi.LocalIndex()];
                    auto found = m_fine_map.find(gi);
                    if (found != m_fine_map.end()) {
                        found->second[idim] = &(m_fine_fluxes[idim][mfi]);
                    } else {
                        Array<FArrayBox*,AMREX_SPACEDIM> t{AMREX_D_DECL(nullptr,nullptr,nullptr)};
                        t[idim] = &(m_fine_fluxes[idim][mfi]);
                        m_fine_map.insert(std::make_pair(gi,t));
                    }
                }
            }
        }
    }

    // For a coarse Box, there are at most two faces per direction abutting fine level.
    {
        BoxArray const fba_coarsened = amrex::coarsen(fba,ref_ratio);
        Array<BoxList,AMREX_SPACEDIM> bl;
        Array<Vector<int>,AMREX_SPACEDIM> procmap;
        Array<Vector<int>,AMREX_SPACEDIM> my_global_indices;
        std::vector<std::pair<int,Box> > isects;
        const std::vector<IntVect>& pshifts = m_crse_geom.periodicity().shiftIntVect();
        Vector<std::pair<Orientation,Box> > cell_faces;
        Vector<Orientation> crsefine_faces;
        for (int i = 0, ncrses = cba.size(); i < ncrses; ++i)
        {
            Box const& ccbx = cba[i];
            Box const& ccbxg1 = amrex::grow(ccbx,1);
            cell_faces.clear();
            crsefine_faces.clear();
            for (OrientationIter orit; orit.isValid(); ++orit) {
                const Orientation ori = orit();
                cell_faces.emplace_back(std::make_pair(ori,amrex::adjCell(ccbx,ori)));
            }
            for (auto const& shift : pshifts) {
                if (!cell_faces.empty()) {
                    fba_coarsened.intersections(ccbxg1+shift, isects);
                    for (auto const& isect : isects) {
                        Box const& b = isect.second-shift;
                        if (ccbx.intersects(b)) {
                            // This coarse box is covered by fine
                            cell_faces.clear();
                            crsefine_faces.clear();
                        } else {
                            auto found = std::find_if(cell_faces.begin(), cell_faces.end(),
                                                      [&] (std::pair<Orientation,Box> const& x)
                                                      { return x.second.contains(b); });
                            if (found != cell_faces.end()) {
                                crsefine_faces.push_back(found->first);
                                cell_faces.erase(found);
                            }
                        }
                    }
                }
            }
            for (auto const& face : crsefine_faces) {
                const int dir = face.coordDir();
                bl[dir].push_back(amrex::bdryNode(ccbx,face));
                procmap[dir].push_back(cdm[i]);
                if (cdm[i] == myproc) my_global_indices[dir].push_back(i);
            }
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (bl[idim].isNotEmpty()) {
                BoxArray faces_ba(std::move(bl[idim]));
                DistributionMapping faces_dm(std::move(procmap[idim]));
                m_crse_fluxes[idim].define(faces_ba, faces_dm, m_ncomp, 0);
                for (MFIter mfi(m_crse_fluxes[idim]); mfi.isValid(); ++mfi) {
                    const int gi = my_global_indices[idim][mfi.LocalIndex()];
                    const Orientation::Side side = (mfi.validbox() == amrex::bdryLo(cba[gi],idim))
                        ? Orientation::low : Orientation::high;
                    const int index = Orientation(idim,side);
                    auto found = m_crse_map.find(gi);
                    if (found != m_crse_map.end()) {
                        found->second[index] = &(m_crse_fluxes[idim][mfi]);
                    } else {
                        Array<FArrayBox*,2*AMREX_SPACEDIM> t{AMREX_D_DECL(nullptr,nullptr,nullptr),
                                                             AMREX_D_DECL(nullptr,nullptr,nullptr)};
                        t[index] = &(m_crse_fluxes[idim][mfi]);
                        m_crse_map.insert(std::make_pair(gi,t));
                    }
                }
            }
        }
    }
}

void FlashFluxRegister::store (int fine_global_index, int dir, FArrayBox const& fine_flux, Real sf)
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_fine_map.find(fine_global_index);
    if (found != m_fine_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,AMREX_SPACEDIM> const& fab_a = found->second;
        if (fab_a[dir]) {
            Box const& b = fab_a[dir]->box();
            Array4<Real> const& dest = fab_a[dir]->array();
            Array4<Real const> const& src = fine_flux.const_array();
            if (dir == 0) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
#if (AMREX_SPACEDIM == 1)
                    dest(i,0,0,n) = src(2*i,0,0,n)*sf;
#elif (AMREX_SPACEDIM == 2)
                    dest(i,j,0,n) = (src(2*i,2*j  ,0,n) +
                                     src(2*i,2*j+1,0,n)) * (0.5*sf);
#elif (AMREX_SPACEDIM == 3)
                    dest(i,j,k,n) = (src(2*i,2*j  ,2*k  ,n) +
                                     src(2*i,2*j+1,2*k  ,n) +
                                     src(2*i,2*j  ,2*k+1,n) +
                                     src(2*i,2*j+1,2*k+1,n)) * (0.25*sf);
#endif
                });
            } else if (dir == 1) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
#if (AMREX_SPACEDIM == 2)
                    dest(i,j,0,n) = (src(2*i  ,2*j,0,n) +
                                     src(2*i+1,2*j,0,n)) * (0.5*sf);
#elif (AMREX_SPACEDIM == 3)
                    dest(i,j,k,n) = (src(2*i  ,2*j,2*k  ,n) +
                                     src(2*i+1,2*j,2*k  ,n) +
                                     src(2*i  ,2*j,2*k+1,n) +
                                     src(2*i+1,2*j,2*k+1,n)) * (0.25*sf);
#endif
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
#if (AMREX_SPACEDIM == 3)
                    dest(i,j,k,n) = (src(2*i  ,2*j  ,2*k,n) +
                                     src(2*i+1,2*j  ,2*k,n) +
                                     src(2*i  ,2*j+1,2*k,n) +
                                     src(2*i+1,2*j+1,2*k,n)) * (0.25*sf);
#endif
                });
            }
        }
    }
}

void FlashFluxRegister::communicate ()
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_crse_fluxes[idim].ParallelCopy(m_fine_fluxes[idim], m_crse_geom.periodicity());
    }
}

void FlashFluxRegister::load (int crse_global_index, int dir, FArrayBox& crse_flux, Real sf) const
{
    AMREX_ASSERT(dir < AMREX_SPACEDIM);
    auto found = m_crse_map.find(crse_global_index);
    if (found != m_crse_map.end()) {
        const int ncomp = m_ncomp;
        Array<FArrayBox*,2*AMREX_SPACEDIM> const& fab_a = found->second;
        Array4<Real> const& dest = crse_flux.array();
        for (int index = dir; index < 2*AMREX_SPACEDIM; index += AMREX_SPACEDIM) {
            if (fab_a[index]) {
                Box const& b = fab_a[index]->box();
                Array4<Real const> const& src = fab_a[index]->const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (b, ncomp, i, j, k, n,
                {
                    dest(i,j,k,n) = src(i,j,k,n) * sf;
                });
            }
        }
    }
}

}
