



int readParms(Parms *parms, int iproc) {

    int ret;

    FILE *fp;
    if(iproc==0)  fp = fopen("./ExecutionParms","w");
    
    ret = getParm("RUN_MODE", &(parms->RUN_MODE), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'RUN_MODE', assuming it is RUN_TWI_REG. \n\n"); parms->RUN_MODE=1; }
    if(parms->RUN_MODE==10)
    {
        ret = getParm("OF_vel_ref_step", &(parms->OF_vel_ref_step), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n Cannot read parameter 'OF_vel_ref_step'. \n\n"); return ret; }
        ret = getParm("OF_vel_mis_step", &(parms->OF_vel_mis_step), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n Cannot read parameter 'OF_vel_mis_step'. \n\n"); return ret; }

        ret = getParm("OF_vel_ref_mag", &(parms->OF_vel_ref_mag), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n Cannot read parameter 'OF_vel_ref_mag'. \n\n"); return ret; }
        ret = getParm("OF_vel_mis_mag", &(parms->OF_vel_mis_mag), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n Cannot read parameter 'OF_vel_mis_mag'. \n\n"); return ret; }
    }

    ret = getParm("computeOnlyGrad_NoLineSearch", &(parms->computeOnlyGrad_NoLineSearch), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'computeOnlyGrad_NoLineSearch', assuming it is zero (do the whole thing). \n\n"); parms->computeOnlyGrad_NoLineSearch=0; }
    

    ret = getStringParm(inputFileName, "executionDirectory", (parms->executionDirectory), "mandatory");
    if(ret==1)  { printf("\n\n Cannot read parameter 'executionDirectory'. \n\n"); return ret; }

    ret = getParm("LX", &(parms->LX), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'LX'. \n\n"); return ret; }

    ret = getParm("LY", &(parms->LY), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'LY'. \n\n"); return ret; }
    
    ret = getParm("LZ", &(parms->LZ), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'LZ'. \n\n"); return ret; }

    ret = getParm("LT", &(parms->LT), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'LT'. \n\n"); return ret; }

    ret = getParm("dx", &(parms->dx), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'dx'. \n\n"); return ret; }

    ret = getParm("dz", &(parms->dz), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'dz'. \n\n"); return ret; }

    ret = getParm("dt", &(parms->dt), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'dt'. \n\n"); return ret; }

    ret = getParm("ox", &(parms->ox), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'ox'. \n\n"); return ret; }

    ret = getParm("oz", &(parms->oz), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'oz'. \n\n"); return ret; }

    ret = getParm("jt_mig", &(parms->jt), "mandatory", TYPE_INT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'jt_mig'. \n\n"); return ret; }

    if(parms->RUN_MODE!=RUN_MODELING) // Perform modeling
    {    
        ret = getParm("EIC", &(parms->jt_EIC), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIC'. Setting it to zero. \n\n");  parms->jt_EIC=0;}

        // Focusing operator and inversion
        ret = getParm("n_iter", &(parms->n_iter), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'n_iter'. \n\n"); return ret; }

        ret = getParm("alphaTry", &(parms->alpha), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'alphaTry'. \n\n"); return ret; }

        ret = getParm("muFocus", &(parms->muFocus), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'muFocus'. \n\n"); return ret; }
    
        ret = getParm("epsFocus", &(parms->epsFocus), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'epsFocus'. \n\n"); return ret; }

        ret = getParm("lambdaNull", &(parms->lambdaNull), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'lambdaNull', assuming it is zero. \n\n"); parms->lambdaNull=0; }

        ret = getParm("tauNull", &(parms->EIMSH_tauNull), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'tauNull', assuming it is zero. \n\n"); parms->EIMSH_tauNull=0; }


        // CIG Filtering
        ret = getParm("applyBandass2EImage", &(parms->applyBandass2EImage), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'applyBandass2EImage', assuming it is zero (not apply bandass to extended image). \n\n"); parms->applyBandass2EImage=0; }
        if(parms->applyBandass2EImage==1)
        { 
            ret = getParm("lambda1", &(parms->lambda1), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'lambda1'. \n\n"); return ret; }
            ret = getParm("lambda2", &(parms->lambda2), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'lambda2'. \n\n"); return ret; }
            ret = getParm("lambda3", &(parms->lambda3), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'lambda3'. \n\n"); return ret; }
            ret = getParm("lambda4", &(parms->lambda4), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'lambda4'. \n\n"); return ret; }
        }
        ret = getParm("eimageTaper_applyXTaper", &(parms->eimageTaper_applyXTaper), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'eimageTaper_applyXTaper', assuming it is zero (not apply X Taper). \n\n"); parms->eimageTaper_applyXTaper=0; }
        if(parms->eimageTaper_applyXTaper==1)
        {       
            ret = getParm("eimageTaper_x1", &(parms->eimageTaper_x1), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_x1'. \n\n"); return ret; }
            ret = getParm("eimageTaper_x2", &(parms->eimageTaper_x2), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_x2'. \n\n"); return ret; }
            ret = getParm("eimageTaper_x3", &(parms->eimageTaper_x3), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_x3'. \n\n"); return ret; }
            ret = getParm("eimageTaper_x4", &(parms->eimageTaper_x4), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_x4'. \n\n"); return ret; }
        }
        ret = getParm("eimageTaper_applyLXTaper", &(parms->eimageTaper_applyLXTaper), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'eimageTaper_applyLXTaper', assuming it is zero (not apply LX Taper). \n\n"); parms->eimageTaper_applyLXTaper=0; }
        if(parms->eimageTaper_applyLXTaper==1)
        {
            ret = getParm("eimageTaper_z_i", &(parms->eimageTaper_z_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z_i'. \n\n"); return ret; }
        
            ret = getParm("eimageTaper_z_f", &(parms->eimageTaper_z_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx1_i", &(parms->eimageTaper_lx1_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx1_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx2_i", &(parms->eimageTaper_lx2_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx2_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx3_i", &(parms->eimageTaper_lx3_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx3_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx4_i", &(parms->eimageTaper_lx4_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx4_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx1_f", &(parms->eimageTaper_lx1_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx1_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx2_f", &(parms->eimageTaper_lx2_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx2_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx3_f", &(parms->eimageTaper_lx3_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx3_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_lx4_f", &(parms->eimageTaper_lx4_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_lx4_f'. \n\n"); return ret; }
    
        }

        ret = getParm("eimageTaper_applyTauTaper", &(parms->eimageTaper_applyTauTaper), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'eimageTaper_applyTauTaper', assuming it is zero (not apply LX Taper). \n\n"); parms->eimageTaper_applyTauTaper=0; }
        if(parms->eimageTaper_applyTauTaper==1)
        {
            ret = getParm("eimageTaper_z_i", &(parms->eimageTaper_z_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z_i'. \n\n"); return ret; }
        
            ret = getParm("eimageTaper_z_f", &(parms->eimageTaper_z_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau1_i", &(parms->eimageTaper_tau1_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau1_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau2_i", &(parms->eimageTaper_tau2_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau2_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau3_i", &(parms->eimageTaper_tau3_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau3_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau4_i", &(parms->eimageTaper_tau4_i), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau4_i'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau1_f", &(parms->eimageTaper_tau1_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau1_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau2_f", &(parms->eimageTaper_tau2_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau2_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau3_f", &(parms->eimageTaper_tau3_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau3_f'. \n\n"); return ret; }

            ret = getParm("eimageTaper_tau4_f", &(parms->eimageTaper_tau4_f), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_tau4_f'. \n\n"); return ret; }
        }

        ret = getParm("eimageTaper_applyZTaper", &(parms->eimageTaper_applyZTaper), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'eimageTaper_applyZTaper', assuming it is zero (not apply Z Taper). \n\n"); parms->eimageTaper_applyZTaper=0; }
        if(parms->eimageTaper_applyZTaper==1)
        {
            ret = getParm("eimageTaper_z0", &(parms->eimageTaper_z0), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z0'. \n\n"); return ret; }
        
            ret = getParm("eimageTaper_z1", &(parms->eimageTaper_z1), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z1'. \n\n"); return ret; }

            ret = getParm("eimageTaper_z2", &(parms->eimageTaper_z2), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z2'. \n\n"); return ret; }

            ret = getParm("eimageTaper_z3", &(parms->eimageTaper_z3), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_z3'. \n\n"); return ret; }
        }

        ret = getParm("eimageFocus_zMin", &(parms->eimageFocus_zMin), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_zMin'. \n\n"); return ret; }
    
        ret = getParm("eimageFocus_zMax", &(parms->eimageFocus_zMax), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_zMax'. \n\n"); return ret; }

        ret = getParm("eimageTaper_applyTaperAngle", &(parms->eimageTaper_applyTaperAngle), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'eimageTaper_applyTaperAngler', \
                               assuming it is zero (not apply angle taper). \n\n"); parms->eimageTaper_applyTaperAngle=0; }
        if(parms->eimageTaper_applyTaperAngle==1)
        {        
            ret = getParm("eimageTaper_theta1", &(parms->eimageTaper_theta1), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_theta1'. \n\n"); return ret; }

            ret = getParm("eimageTaper_theta2", &(parms->eimageTaper_theta2), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_theta2'. \n\n"); return ret; }

            ret = getParm("eimageTaper_theta3", &(parms->eimageTaper_theta3), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_theta3'. \n\n"); return ret; }

            ret = getParm("eimageTaper_theta4", &(parms->eimageTaper_theta4), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'eimageTaper_theta4'. \n\n"); return ret; }
        }

        // Gradient Z tapering
        ret = getParm("gradProc_applyBandpass", &(parms->gradProc_applyBandpass), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'gradProc_applyBandpass', assuming it is zero (not apply bandpass to gradient). \n\n"); parms->gradProc_applyBandpass=0; }
        if(parms->gradProc_applyBandpass==1)
        {        
            ret = getParm("gradProc_LBDA1", &(parms->gradProc_LBDA1), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_LBDA1'. \n\n"); return ret; }
            ret = getParm("gradProc_LBDA2", &(parms->gradProc_LBDA2), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_LBDA2'. \n\n"); return ret; }
            ret = getParm("gradProc_LBDA3", &(parms->gradProc_LBDA3), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_LBDA3'. \n\n"); return ret; }
            ret = getParm("gradProc_LBDA4", &(parms->gradProc_LBDA4), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_LBDA4'. \n\n"); return ret; }
        }


        // Gradient Z tapering
        ret = getParm("gradProc_applyZTaper", &(parms->gradProc_applyZTaper), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'gradProc_applyZTaper', assuming it is zero (not apply Z Taper). \n\n"); parms->gradProc_applyZTaper=0; }
        if(parms->gradProc_applyZTaper==1)
        {        
            ret = getParm("gradProc_z1", &(parms->gradProc_z1), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_z1'. \n\n"); return ret; }

            ret = getParm("gradProc_z2", &(parms->gradProc_z2), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_z2'. \n\n"); return ret; }

            ret = getParm("gradProc_z3", &(parms->gradProc_z3), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_z3'. \n\n"); return ret; }

            ret = getParm("gradProc_z4", &(parms->gradProc_z4), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_z4'. \n\n"); return ret; }
        }
        // Gradient X tapering
        ret = getParm("gradProc_applyXTaper", &(parms->gradProc_applyXTaper), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'gradProc_applyXTaper', assuming it is zero (not apply X Taper). \n\n"); parms->gradProc_applyXTaper=0; }
        if(parms->gradProc_applyXTaper==1)
        {        
            ret = getParm("gradProc_x1", &(parms->gradProc_x1), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_x1'. \n\n"); return ret; }

            ret = getParm("gradProc_x2", &(parms->gradProc_x2), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_x2'. \n\n"); return ret; }

            ret = getParm("gradProc_x3", &(parms->gradProc_x3), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_x3'. \n\n"); return ret; }

            ret = getParm("gradProc_x4", &(parms->gradProc_x4), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_x4'. \n\n"); return ret; }
        }
        // Gradient smoothing
        ret = getParm("gradProc_applySmoothing", &(parms->gradProc_applySmoothing), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'gradProc_applySmoothing', assuming it is zero (not apply smoothing). \n\n"); parms->gradProc_applySmoothing=0; }
        if(parms->gradProc_applySmoothing==1)
        {
            ret = getParm("gradProc_halfLengthZ", &(parms->gradProc_halfLengthZ), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_halfLengthZ'. \n\n"); return ret; }

            ret = getParm("gradProc_halfLengthX", &(parms->gradProc_halfLengthX), "mandatory", TYPE_FLOAT);
            if(ret!=0)  { printf("\n\n ERROR: Cannot read parameter 'gradProc_halfLengthX'. \n\n"); return ret; }
        }
        
        // Angle of dip
        ret = getParm("EIMSH_dipMin", &(parms->EIMSH_dipMin), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_dipMin', assuming it is zero. \n\n"); parms->EIMSH_dipMin = 0.0; }
        ret = getParm("EIMSH_dipMax", &(parms->EIMSH_dipMax), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_dipMax', assuming it is zero. \n\n"); parms->EIMSH_dipMax = 0.0; }
        ret = getParm("EIMSH_ddip", &(parms->EIMSH_ddip), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_ddip', assuming it is 1.0. \n\n");    parms->EIMSH_ddip   = 1.0; }
        
        ret = getParm("EIMSH_lambdaMax", &(parms->EIMSH_lambdaMax), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'EIMSH_lambdaMax'. \n\n"); return ret; }
        ret = getParm("EIMSH_dlambda", &(parms->EIMSH_dlambda), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'EIMSH_dlambda'. \n\n"); return ret; }

        ret = getParm("EIMSH_CIGS_TAU", &(parms->EIMSH_CIGS_TAU), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_CIGS_TAU', assuming it is 0. \n\n"); parms->EIMSH_CIGS_TAU = 0; }
        if(parms->EIMSH_CIGS_TAU) 
        {
            ret = getParm("EIMSH_tauMax", &(parms->EIMSH_tauMax), "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'EIMSH_tauMax'. \n\n"); return ret; }
            ret = getParm("EIMSH_dtau"  , &(parms->EIMSH_dtau)  , "mandatory", TYPE_FLOAT);
            if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'EIMSH_dtau'.   \n\n"); return ret; }
        }

        ret = getParm("EIMSH_CIGS_ORTH", &(parms->EIMSH_CIGS_ORTH), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_CIGS_ORTH', assuming it is 1. \n\n"); parms->EIMSH_CIGS_ORTH = 1; }

        if(parms->EIMSH_CIGS_ORTH) {
            ret = getParm("EIMSH_CIGS_VCIG", &(parms->EIMSH_CIGS_VCIG), "optional", TYPE_INT);
            if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_CIGS_VCIG', assuming it is 0. \n\n"); parms->EIMSH_CIGS_VCIG = 0; }
        }

        ret = getParm("EIMSH_focusingQC", &(parms->EIMSH_focusingQC), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_focusingQC', assuming it is 0. \n\n"); parms->EIMSH_focusingQC = 0; }


        ret = getParm("EIMSH_thetaMax", &(parms->EIMSH_thetaMax), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_thetaMax', assmuming EIMSH_thetaMax = 45.0'. \n\n"); parms->EIMSH_thetaMax=45.0; }

        ret = getParm("EIMSH_dtheta", &(parms->EIMSH_dtheta), "optional", TYPE_FLOAT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'EIMSH_dtheta', assmuming EIMSH_dtheta = 2.0'. \n\n");  parms->EIMSH_thetaMax=2.0; }

    }

    ret = getParm("f1", &(parms->f1), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'f1'. \n\n"); return ret; }

    ret = getParm("f2", &(parms->f2), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'f2'. \n\n"); return ret; }

    ret = getParm("f3", &(parms->f3), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'f3'. \n\n"); return ret; }

    ret = getParm("f4", &(parms->f4), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'f4'. \n\n"); return ret; }

    ret = getParm("offMax", &(parms->offMax), "optional", TYPE_FLOAT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'offMax', assuming it is infinity. \n\n"); parms->offMax=1.0e12; }

    ret = getParm("offPerc", &(parms->offPerc), "offPerc", TYPE_FLOAT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'offPerc', assuming it is zero. \n\n"); parms->offPerc=0.0; }
    
    ret = getParm("acquisitionGeom", &(parms->acquisitionGeom), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'acquisitionGeom', assuming it is zero (split-spread). \n\n"); parms->acquisitionGeom=0; }
    
    
    ret = getParm("ngpu", &(parms->ngpu), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'ngpu', assuming it is 1. \n\n"); parms->ngpu=1; }
    ret = getParm("firstGPU", &(parms->firstGPU), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'firstGPU', assuming it is 0. \n\n"); parms->firstGPU=0; }

    
    ret = getParm("nxInput_vel", &(parms->nxInput_vel), "mandatory", TYPE_INT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'nxInput_vel'. \n\n"); return ret; }
    ret = getParm("nzInput_vel", &(parms->nzInput_vel), "mandatory", TYPE_INT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'nzInput_vel'. \n\n"); return ret; }
    ret = getParm("deltaXInput_vel", &(parms->deltaXInput_vel), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'deltaXInput_vel'. \n\n"); return ret; }
    ret = getParm("deltaZInput_vel", &(parms->deltaZInput_vel), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'deltaZInput_vel'. \n\n"); return ret; }
    ret = getParm("OrigXInput_vel", &(parms->OrigXInput_vel), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'OrigXInput_vel'. \n\n"); return ret; }
    ret = getParm("OrigZInput_vel", &(parms->OrigZInput_vel), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'OrigZInput_vel'. \n\n"); return ret; }
    ret = getParm("velScalar", &(parms->velScalar), "mandatory", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'velScalar'. \n\n"); return ret; }

    ret = getStringParm(inputFileName, "filePathVel", (parms->filePathVel), "mandatory");
    if(ret==1)  { printf("\n\n Cannot read parameter 'filePathVel'. \n\n"); return ret; }

    if(parms->RUN_MODE==RUN_TWI_REG || parms->RUN_MODE==RUN_TWI_OBJFUNCEXP) // Perform regular TWI or Objective Fucntion Evaluation
    {
        ret = getStringParm(inputFileName, "inputSeismicFileName", (parms->inputSeismicFileName), "mandatory");
        if(ret==1)  { printf("\n\n Cannot read parameter 'inputSeismicFileName'. \n\n"); return ret; }
    }
    if(parms->RUN_MODE==RUN_TWI_TL) // Perform TL-TWI
    {
        ret = getStringParm(inputFileName, "inputMoniSeismicFileName", (parms->inputSeismicFileName_moni), "mandatory");
        if(ret==1)  { printf("\n\n Cannot read parameter 'inputMoniSeismicFileName'. \n\n"); return ret; }

        ret = getStringParm(inputFileName, "inputBaseSeismicFileName", (parms->inputSeismicFileName_base), "mandatory");
        if(ret==1)  { printf("\n\n Cannot read parameter 'inputBaseSeismicFileName'. \n\n"); return ret; }
    }

    ret = getParm("performEBM", &(parms->performEBM), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'performEBM', assuming it is zero. \n\n"); parms->performEBM=0; }


    // ==============================================================
    // ======================  New parameters  ======================
    // ==============================================================
    
    // Parameters for top mute of shot
    ret = getParm("TopMute_apply", &(parms->TopMute_apply), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'TopMute_apply', assuming it is zero (do not apply top mute). \n\n"); parms->TopMute_apply=0; }
    if(parms->TopMute_apply)
    {
        ret = getParm("TopMute_t0", &(parms->TopMute_t0), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'TopMute_t0'. \n\n"); return ret; }
        ret = getParm("TopMute_t1", &(parms->TopMute_t1), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'TopMute_t1'. \n\n"); return ret; }

        ret = getParm("TopMute_off0", &(parms->TopMute_off0), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'TopMute_off0'. \n\n"); return ret; }
        ret = getParm("TopMute_off1", &(parms->TopMute_off1), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'TopMute_off1'. \n\n"); return ret; }

        ret = getParm("TopMute_ramp", &(parms->TopMute_ramp), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'TopMute_ramp'. \n\n"); return ret; }
    }

    // Parameters for bottom mute of shot
    ret = getParm("BottomMute_apply", &(parms->BottomMute_apply), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'BottomMute_apply', assuming it is zero (do not apply top mute). \n\n"); parms->BottomMute_apply=0; }
    if(parms->BottomMute_apply)
    {
        ret = getParm("BottomMute_t0", &(parms->BottomMute_t0), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'BottomMute_t0'. \n\n"); return ret; }
        ret = getParm("BottomMute_t1", &(parms->BottomMute_t1), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'BottomMute_t1'. \n\n"); return ret; }

        ret = getParm("BottomMute_off0", &(parms->BottomMute_off0), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'BottomMute_off0'. \n\n"); return ret; }
        ret = getParm("BottomMute_off1", &(parms->BottomMute_off1), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'BottomMute_off1'. \n\n"); return ret; }

        ret = getParm("BottomMute_ramp", &(parms->BottomMute_ramp), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'BottomMute_ramp'. \n\n"); return ret; }
    }

    ret = getParm("normalizeInputDataAmplitudes", &(parms->normalizeInputDataAmplitudes), "optional", TYPE_INT);
    if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'normalizeInputDataAmplitudes', assuming not normalize amplitudes of input seismic data. \n\n"); parms->normalizeInputDataAmplitudes=0; }
    if(iproc==0)  fprintf(fp, "\n  Parameter normalizeInputDataAmplitudes:    <%d>", parms->normalizeInputDataAmplitudes);

    ret = getParm("AppertureExtension_Iline", &(parms->AppertureExtension_x), "optional", TYPE_FLOAT);
    if(ret==1)  { printf("\n\n WARNING: Cannot read parameter 'AppertureExtension_Iline', assuming AppertureExtension_Iline=0.0 (not adding X apperture) \n\n"); parms->AppertureExtension_x=0.0; }
    if(iproc==0)  fprintf(fp, "\n  Parameter AppertureExtension_x:    <%f>", parms->AppertureExtension_x);

    ret = getParm("nxb", &(parms->nxb), "optional", TYPE_INT);
    if(ret==1)  { printf("\n\n WARNING: Cannot read parameter 'nxb'. Assuming nxb=40 \n\n"); parms->nxb = 40; }
    if(iproc==0)  fprintf(fp, "\n  Parameter nxb:  < %d >", parms->nxb);

    ret = getParm("nyb", &(parms->nyb), "optional", TYPE_INT);
    if(ret==1)  { printf("\n\n WARNING: Cannot read parameter 'nyb'. Assuming nyb=40 \n\n"); parms->nyb = 40; }
    if(iproc==0)  fprintf(fp, "\n  Parameter nyb:  < %d >", parms->nyb);

    ret = getParm("nzb", &(parms->nzb), "optional", TYPE_INT);
    if(ret==1)  { printf("\n\n WARNING: Cannot read parameter 'nzb'. Assuming nzb=40 \n\n"); parms->nzb = 40; }
    if(iproc==0)  fprintf(fp, "\n  Parameter nzb:  < %d >", parms->nzb);


    if(parms->RUN_MODE==0) // Perform modeling
    {   
        ret = getStringParm(inputFileName, "outputSeismicFileName", (parms->outputSeismicFileName), "mandatory");
        if(ret==1)  { printf("\n\n Cannot read parameter 'outputSeismicFileName'. \n\n"); return ret; }

        ret = getStringParm(inputFileName, "filePathDAVel", (parms->filePathDAVel), "mandatory");
        if(ret==1)  { printf("\n\n Cannot read parameter 'filePathDAVel'. \n\n"); return ret; }

        // Shot geom parms
        ret = getParm("shotOX", &(parms->shotOX), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotOX'. \n\n"); return ret; }
    
        ret = getParm("shotOY", &(parms->shotOY), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotOY'. \n\n"); return ret; }

        ret = getParm("shotOZ", &(parms->shotOZ), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotOZ'. \n\n"); return ret; }

        ret = getParm("shotDX", &(parms->shotDX), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotDX'. \n\n"); return ret; }
    
        ret = getParm("shotDY", &(parms->shotDY), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotDY'. \n\n"); return ret; }

        ret = getParm("shotDZ", &(parms->shotDZ), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotDZ'. \n\n"); return ret; }

        ret = getParm("shotNX", &(parms->shotNX), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotNX'. \n\n"); return ret; }
    
        ret = getParm("shotNY", &(parms->shotNY), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotNY'. \n\n"); return ret; }

        ret = getParm("shotNZ", &(parms->shotNZ), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'shotNZ'. \n\n"); return ret; }
    
        // Rec geom parms
        ret = getParm("recOX", &(parms->recOX), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recOX'. \n\n"); return ret; }
    
        ret = getParm("recOY", &(parms->recOY), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recOY'. \n\n"); return ret; }

        ret = getParm("recOZ", &(parms->recOZ), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recOZ'. \n\n"); return ret; }

        ret = getParm("recDX", &(parms->recDX), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recDX'. \n\n"); return ret; }
    
        ret = getParm("recDY", &(parms->recDY), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recDY'. \n\n"); return ret; }

        ret = getParm("recDZ", &(parms->recDZ), "mandatory", TYPE_FLOAT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recDZ'. \n\n"); return ret; }

        ret = getParm("recNX", &(parms->recNX), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recNX'. \n\n"); return ret; }
    
        ret = getParm("recNY", &(parms->recNY), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recNY'. \n\n"); return ret; }

        ret = getParm("recNZ", &(parms->recNZ), "mandatory", TYPE_INT);
        if(ret==1)  { printf("\n\n ERROR: Cannot read parameter 'recNZ'. \n\n"); return ret; }

        ret = getParm("recCoordRefFromShot", &(parms->recCoordRefFromShot), "optional", TYPE_INT);
        if(ret!=0)  { printf("\n\n WARNING: Cannot read parameter 'recCoordRefFromShot'. Assuming recCoordRefFromShot=1 \n\n"); parms->recCoordRefFromShot = 1; }
        if(iproc==0)  fprintf(fp, "\n  Parameter recCoordRefFromShot:  < %d >", parms->recCoordRefFromShot);
    }

return 0;
}

int getParm(char *parmName, void* parm, char *requirement, parmType TYPE)
{
    int  stat=0;
    int  read=0;
    char inputParm[1024];

    FILE *fp;
    fp = fopen(inputFileName,"r");
    if(fp==NULL)
    {
        printf("\n\n Cannot open file with input parameters. \
                     file name: %s \n\n", inputFileName);
        stat = 2;
        return stat;
    }
    while(!feof(fp))
    {
        fscanf(fp, "%s", inputParm);
        if( (strncmp(inputParm, parmName, strlen(parmName))) == 0)
        {
            if(TYPE==TYPE_INT)
                fscanf(fp,"%d" , (int    *) parm);
            if(TYPE==TYPE_FLOAT)
                fscanf(fp,"%f" , (float  *) parm);
            if(TYPE==TYPE_DOUBLE)
                fscanf(fp,"%lf", (double *) parm);
            if(TYPE==TYPE_STRING) {
                // printf("\n\n  >>>>>>>>>>>>>>>>>>>> !!!!!!!!!!!  TYPE_STRING !!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<< \n\n");
                fscanf(fp,"%s" , (char   *) parm);
                // printf("\n\n  parms=%s \n\n"  , parm);
            }
                
            read++;
        }
    }
    fclose(fp);

    if(read < 0)
        stat = -1;
    if(read == 0)
    {
        if( strcmp(requirement,"mandatory") == 0 )
            stat = 1;
        if( strcmp(requirement,"optional") == 0 )
            stat = -1;
    }
    if(read >= 1)
        stat = 0;

    
return stat;
}

int getStringParm(char *inputFileName, char *parmName, char* parm, char *requirement)
{
    int  stat=0;
    int  read=0;
    char inputParm[1024];
    char inputParmValue[1024];

    FILE *fp;
    fp = fopen(inputFileName,"r");
    if(fp==NULL)
    {
        printf("\n\n Cannot open file with input parameters. \
                     file name: %s \n\n", inputFileName);
        stat = 2;
        return stat;
    }
    while(!feof(fp))
    {
        fscanf(fp, "%s", inputParm);
        if( (strncmp(inputParm, parmName, strlen(parmName))) == 0)
        {
            fscanf(fp, "%s", inputParmValue);
            read++;
        }
    }
    fclose(fp);

    if(read < 0)
        stat = -1;
    if(read == 0)
    {
        if( strcmp(requirement,"mandatory") == 0 )
            stat = 1;
        if( strcmp(requirement,"optional") == 0 )
            stat = -1;
    }
    if(read >= 1)
    { 
        strncpy(parm, inputParmValue, 1024);
        stat = 0;
    }
        

    
return stat;
}

int readCommandLineArg(void* parm, int argc, char** argv, char *argName, char *requirement, parmType TYPE)
{ 
    int stat, read=0;

    for (int i = 0; i < argc; ++i) 
    {

        // Initialize arrays of characters
        char inputParm[1024];
        char inputValue[1024];
        memset(inputParm,  0, sizeof inputParm);
        memset(inputValue, 0, sizeof inputValue);

        // Reads argument name
        int ic=0;
        while(argv[i][ic]!='\0')
        {
            if(argv[i][ic]=='=') { 
                ic++;
                break;
            }
            else {
                inputParm[ic] = argv[i][ic];
                ic++;
            }
        }

        // Checks if this is the requested argument
        int nc = ic;
        if(strncmp(argName, inputParm, ic-1)==0)
        {
            // printf("\n\n  found parm %s   inputParm=%s\n\n", argName, inputParm);
            while(argv[i][ic]!='\0')
            {
                inputValue[ic-nc] = argv[i][ic];
                ic++;
            }
            if(TYPE==TYPE_INT)        sscanf(inputValue, "%d" , (int    *) parm);
            if(TYPE==TYPE_FLOAT)      sscanf(inputValue, "%f" , (float  *) parm);
            if(TYPE==TYPE_DOUBLE)     sscanf(inputValue, "%lf", (double *) parm);
            if(TYPE==TYPE_STRING)     sscanf(inputValue, "%s" , (char   *) parm);
            read = 1;
        }
    }

    if(read == 1) 
    {
        stat = 0;
    }    
    if(read == 0)
    {
        if( strcmp(requirement,"mandatory") == 0 )    
        stat = 1;
        
        if( strcmp(requirement,"optional" ) == 0 )    
        stat = 0;
    }

return stat;
}