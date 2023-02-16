#include "write_proj_band_lcao.h"
#include "write_orb_info.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"

void ModuleIO::write_proj_band_lcao(const psi::Psi<double> *psid,
    const psi::Psi<std::complex<double>> *psi,
    LCAO_Hamilt& uhm,
    const elecstate::ElecState* pelec,
    const K_Vectors* kv,
    const UnitCell &ucell,
    const LCAO_Orbitals &ORB,
    Grid_Driver &GridD)
{
    ModuleBase::TITLE("ModuleIO","write_proj_band_lcao");
    ModuleBase::timer::tick("ModuleIO", "write_proj_band_lcao");

    const Parallel_Orbitals* pv = uhm.LM->ParaV;

    int nspin0 = 1;
    if (GlobalV::NSPIN == 2)
    nspin0 = 2;
    int nks = 0;
    if (nspin0 == 1)
    {
        nks = kv->nkstot;
    }
    else if (nspin0 == 2)
    {
        nks = kv->nkstot / 2;
    }

    ModuleBase::ComplexMatrix weightk;
    ModuleBase::matrix weight;
    int NUM = 0;
    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        NUM = GlobalV::NLOCAL * GlobalV::NBANDS * nspin0;
        weightk.create(nspin0, GlobalV::NBANDS * GlobalV::NLOCAL, true);
        weight.create(nspin0, GlobalV::NBANDS * GlobalV::NLOCAL, true);
    }
    else
    {
        NUM = GlobalV::NLOCAL * GlobalV::NBANDS * kv->nks;
        weightk.create(kv->nks, GlobalV::NBANDS * GlobalV::NLOCAL, true);
        weight.create(kv->nks, GlobalV::NBANDS * GlobalV::NLOCAL, true);
    }

    for (int is = 0; is < nspin0; is++)
    {

        if (GlobalV::GAMMA_ONLY_LOCAL)
        {
            std::vector<ModuleBase::matrix> Mulk;
            Mulk.resize(1);
            Mulk[0].create(pv->ncol, pv->nrow);

            psid->fix_k(is);
            for (int i = 0; i < GlobalV::NBANDS; ++i)
            {
                const int NB = i + 1;

                const double one_float = 1.0, zero_float = 0.0;
                const int one_int = 1;
#ifdef __MPI
                const char T_char = 'T';
                pdgemv_(&T_char,
                        &GlobalV::NLOCAL,
                        &GlobalV::NLOCAL,
                        &one_float,
                        uhm.LM->Sloc.data(),
                        &one_int,
                        &one_int,
                        pv->desc,
                        psid->get_pointer(),
                        &one_int,
                        &NB,
                        pv->desc,
                        &one_int,
                        &zero_float,
                        Mulk[0].c,
                        &one_int,
                        &NB,
                        pv->desc,
                        &one_int);
#endif
                for (int j = 0; j < GlobalV::NLOCAL; ++j)
                {

                    if (pv->in_this_processor(j, i))
                    {

                        const int ir = pv->trace_loc_row[j];
                        const int ic = pv->trace_loc_col[i];
                        weightk(is, i * GlobalV::NLOCAL + j) = Mulk[0](ic, ir) * psid[0](ic, ir);
                    }
                }
            } // ib
        } // if
        else
        {
            GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                             GlobalV::OUT_LEVEL,
                                                             ORB.get_rcutmax_Phi(),
                                                             ucell.infoNL.get_rcutmax_Beta(),
                                                             GlobalV::GAMMA_ONLY_LOCAL);

            atom_arrange::search(GlobalV::SEARCH_PBC,
                                 GlobalV::ofs_running,
                                 GridD,
                                 ucell,
                                 GlobalV::SEARCH_RADIUS,
                                 GlobalV::test_atom_input); // qifeng-2019-01-21

            uhm.LM->allocate_HS_R(pv->nnr);
            uhm.LM->zeros_HSR('S');
            uhm.genH.calculate_S_no(uhm.LM->SlocR.data());
            uhm.genH.build_ST_new('S', false, ucell, uhm.LM->SlocR.data());
            std::vector<ModuleBase::ComplexMatrix> Mulk;
            Mulk.resize(1);
            Mulk[0].create(pv->ncol, pv->nrow);

            for (int ik = 0; ik < kv->nks; ik++)
            {

                if (is == kv->isk[ik])
                {
                    uhm.LM->allocate_HS_k(pv->nloc);
                    uhm.LM->zeros_HSk('S');
                    uhm.LM->folding_fixedH(ik);

                    psi->fix_k(ik);
                    psi::Psi<std::complex<double>> Dwfc(psi[0], 1);
                    std::complex<double>* p_dwfc = Dwfc.get_pointer();
                    for (int index = 0; index < Dwfc.size(); ++index)
                    {
                        p_dwfc[index] = conj(p_dwfc[index]);
                    }

                    for (int i = 0; i < GlobalV::NBANDS; ++i)
                    {
                        const int NB = i + 1;

                        const double one_float[2] = {1.0, 0.0}, zero_float[2] = {0.0, 0.0};
                        const int one_int = 1;
                        //   const int two_int=2;
                        const char T_char = 'T'; // N_char='N',U_char='U'
#ifdef __MPI
                        pzgemv_(&T_char,
                                &GlobalV::NLOCAL,
                                &GlobalV::NLOCAL,
                                &one_float[0],
                                uhm.LM->Sloc2.data(),
                                &one_int,
                                &one_int,
                                pv->desc,
                                p_dwfc,
                                &one_int,
                                &NB,
                                pv->desc,
                                &one_int,
                                &zero_float[0],
                                Mulk[0].c,
                                &one_int,
                                &NB,
                                pv->desc,
                                &one_int);
#endif
                        for (int j = 0; j < GlobalV::NLOCAL; ++j)
                        {

                            if (pv->in_this_processor(j, i))
                            {

                                const int ir = pv->trace_loc_row[j];
                                const int ic = pv->trace_loc_col[i];

                                weightk(ik, i * GlobalV::NLOCAL + j) = Mulk[0](ic, ir) * psi[0](ic, ir);
                            }
                        }

                    } // ib

                } // if
            } // ik
#ifdef __MPI
            atom_arrange::delete_vector(GlobalV::ofs_running,
                                        GlobalV::SEARCH_PBC,
                                        GridD,
                                        ucell,
                                        GlobalV::SEARCH_RADIUS,
                                        GlobalV::test_atom_input);
#endif
        } // else
#ifdef __MPI
        MPI_Reduce(weightk.real().c, weight.c, NUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

        if (GlobalV::MY_RANK == 0)
        {
            std::stringstream ps2;
            ps2 << GlobalV::global_out_dir << "PBANDS_" << is + 1;
            GlobalV::ofs_running << "\n Output projected bands in file: " << ps2.str() << std::endl;
            std::ofstream out(ps2.str().c_str());

            out << "<pband>" << std::endl;
            out << "<nspin>" << GlobalV::NSPIN << "</nspin>" << std::endl;
            if (GlobalV::NSPIN == 4)
                out << "<norbitals>" << std::setw(2) << GlobalV::NLOCAL / 2 << "</norbitals>" << std::endl;
            else
                out << "<norbitals>" << std::setw(2) << GlobalV::NLOCAL << "</norbitals>" << std::endl;
            out << "<band_structure nkpoints=\"" << nks << "\" nbands=\"" << GlobalV::NBANDS << "\" units=\"eV\">"
                << std::endl;
            if (GlobalV::GAMMA_ONLY_LOCAL)
            {
                for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                    out << " " << (pelec->ekb(is * nks, ib)) * ModuleBase::Ry_to_eV;
                out << std::endl;
            }
            else
            {
                for (int ik = 0; ik < nks; ik++)
                {
                    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                        out << " " << (pelec->ekb(ik + is * nks, ib)) * ModuleBase::Ry_to_eV;
                    out << std::endl;
                }
            }
            out << "</band_structure>" << std::endl;

            for (int i = 0; i < ucell.nat; i++)
            {
                int a = ucell.iat2ia[i];
                int t = ucell.iat2it[i];
                Atom* atom1 = &ucell.atoms[t];
                const int s0 = ucell.itiaiw2iwt(t, a, 0);
                for (int j = 0; j < atom1->nw; ++j)
                {
                    const int L1 = atom1->iw2l[j];
                    const int N1 = atom1->iw2n[j];
                    const int m1 = atom1->iw2m[j];
                    const int w = ucell.itiaiw2iwt(t, a, j);

                    // out << "</energy_values>" <<std::endl;
                    out << "<orbital" << std::endl;
                    out << std::setw(6) << "index=\"" << std::setw(40) << w + 1 << "\"" << std::endl;
                    out << std::setw(5) << "atom_index=\"" << std::setw(40) << i + 1 << "\"" << std::endl;
                    out << std::setw(8) << "species=\"" << ucell.atoms[t].label << "\"" << std::endl;
                    out << std::setw(2) << "l=\"" << std::setw(40) << L1 << "\"" << std::endl;
                    out << std::setw(2) << "m=\"" << std::setw(40) << m1 << "\"" << std::endl;
                    out << std::setw(2) << "z=\"" << std::setw(40) << N1 + 1 << "\"" << std::endl;
                    out << ">" << std::endl;
                    out << "<data>" << std::endl;
                    if (GlobalV::GAMMA_ONLY_LOCAL)
                    {
                        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                        {
                            if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
                                out << std::setw(13) << weight(is, ib * GlobalV::NLOCAL + w);
                            else if (GlobalV::NSPIN == 4)
                            {
                                int w0 = w - s0;
                                out << std::setw(13)
                                    << weight(is, ib * GlobalV::NLOCAL + s0 + 2 * w0)
                                           + weight(is, ib * GlobalV::NLOCAL + s0 + 2 * w0 + 1);
                            }
                        }
                        out << std::endl;
                    }
                    else
                    {
                        for (int ik = 0; ik < nks; ik++)
                        {
                            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
                            {
                                if (GlobalV::NSPIN == 1)
                                    out << std::setw(13) << weight(ik, ib * GlobalV::NLOCAL + w);
                                else if (GlobalV::NSPIN == 2)
                                    out << std::setw(13) << weight(ik + nks * is, ib * GlobalV::NLOCAL + w);
                                else if (GlobalV::NSPIN == 4)
                                {
                                    int w0 = w - s0;
                                    out << std::setw(13)
                                        << weight(ik, ib * GlobalV::NLOCAL + s0 + 2 * w0)
                                               + weight(ik, ib * GlobalV::NLOCAL + s0 + 2 * w0 + 1);
                                }
                            }
                            out << std::endl;
                        }
                    }
                    out << "</data>" << std::endl;
                    out << "</orbital>" << std::endl;
                } // j
            } // i

            out << "</pband>" << std::endl;
            out.close();
        }
    } // is
    ModuleIO::write_orb_info();

    ModuleBase::timer::tick("ModuleIO", "write_proj_band_lcao");
    return;
}
