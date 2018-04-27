function DoScanAtCPU2()
global VSeq;
global VObj;
global VCtl;
global VMag;
global VCoi;
global VVar;
global VSig;
%    const mwSize *SpinMxDims;
%	float *MzBase, *MyBase, *MxBase, *Mz, *My, *Mx, *RhoBase, *T1Base, *T2Base, *Rho, *T1, *T2;

%/* loop control */
i=int32(0); j=int32(0); s=int32(0); Signali=int32(0); Spini=int32(0); Typei=int32(0); RxCoili=int32(0); TxCoili=int32(0);
Slicei=int32(0);
%mwSize MaxStep, MaxutsStep, MaxrfStep, MaxGzStep, MaxGyStep, MaxGxStep, *SpinNum, *TypeNum, *TxCoilNum, *RxCoilNum, *SignalNum;
flag = zeros(1,6,'double');
% /* IPP or FW buffer */
%     Ipp32f *buffer1, *buffer2, *buffer3, *buffer4, *buffer;
% /* Function status */
%     int ExtCall;
% /* assign pointers */
Gyro            = VObj.Gyro;
MzBase          = VObj.Mz;
MyBase          = VObj.My;
MxBase          = VObj.Mx;
RhoBase         = VObj.Rho;
T1Base          = VObj.T1;
T2Base          = VObj.T2;
SpinNum         = VObj.SpinNum;
TypeNum         = VObj.TypeNum;

dB0Base         = VMag.dB0;
dWRndBase       = VMag.dWRnd;
GzgridBase      = VMag.Gzgrid;
GygridBase      = VMag.Gygrid;
GxgridBase      = VMag.Gxgrid;

TxCoilmgBase    = VCoi.TxCoilmg;
TxCoilpeBase    = VCoi.TxCoilpe;
RxCoilx         = VCoi.RxCoilx;
RxCoily         = VCoi.RxCoily;
TxCoilNum       = VCoi.TxCoilNum;
RxCoilNum       = VCoi.RxCoilNum;
TxCoilDefault   = VCoi.TxCoilDefault;
RxCoilDefault   = VCoi.RxCoilDefault;

CS              = VCtl.CS;
TRNum           = VCtl.TRNum;
RunMode         = VCtl.RunMode;
MaxThreadNum    = VCtl.MaxThreadNum;

utsLine         = VSeq.utsLine;
tsLine          = VSeq.tsLine;
rfAmpLine       = VSeq.rfAmpLine;
rfPhaseLine     = VSeq.rfPhaseLine;
rfFreqLine      = VSeq.rfFreqLine;
rfCoilLine      = VSeq.rfCoilLine;
GzAmpLine       = VSeq.GzAmpLine;
GyAmpLine       = VSeq.GyAmpLine;
GxAmpLine       = VSeq.GxAmpLine;
ADCLine         = VSeq.ADCLine;
ExtLine         = VSeq.ExtLine;
flagsLine       = VSeq.flagsLine;

MaxStep         = length(VSeq.tsLine);
MaxutsStep      = length(VSeq.utsLine);
MaxrfStep       = length(VSeq.rfAmpLine);
MaxGzStep       = length(VSeq.GzAmpLine);
MaxGyStep       = length(VSeq.GyAmpLine);
MaxGxStep       = length(VSeq.GxAmpLine);

t               = VVar.t;
dt              = VVar.dt;
rfAmp           = VVar.rfAmp;
rfPhase         = VVar.rfPhase;
rfFreq          = VVar.rfFreq;
rfCoil          = VVar.rfCoil;
rfRef           = VVar.rfRef;
GzAmp           = VVar.GzAmp;
GyAmp           = VVar.GyAmp;
GxAmp           = VVar.GxAmp;
ADC             = VVar.ADC;
Ext             = VVar.Ext;
KzTmp           = VVar.Kz;
KyTmp           = VVar.Ky;
KxTmp           = VVar.Kx;
utsi            = VVar.utsi;
rfi             = VVar.rfi;
Gzi             = VVar.Gzi;
Gyi             = VVar.Gyi;
Gxi             = VVar.Gxi;
ADCi            = VVar.ADCi;
Exti            = VVar.Exti;
TRCount         = VVar.TRCount;

Sy              = VSig.Sy;
Sx              = VSig.Sx;
Kz              = VSig.Kz;
Ky              = VSig.Ky;
Kx              = VSig.Kx;
MzsBase         = VSig.Mz;
MysBase         = VSig.My;
MxsBase         = VSig.Mx;
Muts            = VSig.Muts;
SignalNum       = VSig.SignalNum;

%  /* get dimensions of spin matrix */
SpinMxDimNum    = length(size(VObj.Mz));%mxGetNumberOfDimensions(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "Mz"));
%     SpinMxDims      = (mwSize*) mxCalloc(SpinMxDimNum, sizeof(mwSize));
SpinMxDims      = size(VObj.Mz);%mxGetDimensions(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "Mz"));
SpinMxNum       = prod(SpinMxDims(1:2));%SpinMxDims[0] * SpinMxDims[1];



if (SpinMxDimNum == 2)
    SpinMxSliceNum = 1;
else
    SpinMxSliceNum = SpinMxDims(3);
end

% /* assign spins to multi-threads */
%     if (SpinMxSliceNum < *MaxThreadNum)
%         ThreadNum = SpinMxSliceNum;
%     else
%         ThreadNum = *MaxThreadNum;  /* Full CPU load */
%     /*  Number of parallel threads that the code will try to use (integer).
%         Set to 1 if you want to use a single core, set >1 for multithreading.
%         This number should be less than the #cpus per machine x #cores/cpu x #thread/core.
%         i.e. a machine which two quad-core cpu's could have upto 8 threads.
%         Note that if there are fewer images than threads...
%         then the code will automatically turn down the number of threads
%         (since the extra ones do nothing and waste resources) */
%
% /* set buffer */
%buffer1 = zeros(SpinMxSliceNum,SpinMxNum,'single');
%buffer2 = zeros(SpinMxSliceNum,SpinMxNum,'single');
%buffer3 = zeros(SpinMxSliceNum,SpinMxNum,'single');
%buffer4 = zeros(SpinMxSliceNum,SpinMxNum,'single');
buffer1 = zeros(1, SpinMxSliceNum*SpinMxNum, 'single');
buffer2 = zeros(1, SpinMxSliceNum*SpinMxNum, 'single');
buffer3 = zeros(1, SpinMxSliceNum*SpinMxNum, 'single');
buffer4 = zeros(1, SpinMxSliceNum*SpinMxNum, 'single');
buffer = zeros(1,1,'single');
%
% /* start simulator execution loop */
fprintf("TR Counts: %d of %d\n", 1, TRNum);
while (i < MaxStep-1)
    % /* check MR sequence pulse flag */
    flag(1:6)=0;
    if (tsLine(i+1) ~= tsLine(i+2))
        flag(1)= flag(1) + flagsLine(1,i+1);
        flag(2)= flag(2) + flagsLine(2,i+1);
        flag(3)= flag(3) + flagsLine(3,i+1);
        flag(4)= flag(4) + flagsLine(4,i+1);
        flag(5)= flag(5) + flagsLine(5,i+1);
        flag(6)= flag(6) + flagsLine(6,i+1);
        i = i + 1;
    else
        flag(1)= flag(1) + flagsLine(1,i+1);
        flag(2)= flag(2) + flagsLine(2,i+1);
        flag(3)= flag(3) + flagsLine(3,i+1);
        flag(4)= flag(4) + flagsLine(4,i+1);
        flag(5)= flag(5) + flagsLine(5,i+1);
        flag(6)= flag(6) + flagsLine(6,i+1);
        
        while (tsLine(i+1)==tsLine(i+2))
            flag(1)= flag(1) + flagsLine(1,i+2);
            flag(2)= flag(2) + flagsLine(2,i+2);
            flag(3)= flag(3) + flagsLine(3,i+2);
            flag(4)= flag(4) + flagsLine(4,i+2);
            flag(5)= flag(5) + flagsLine(5,i+2);
            flag(6)= flag(6) + flagsLine(6,i+2);
            i = i + 1;
            if (i==MaxStep-1)
                break;
            end
        end
        i = i + 1;
    end
    
    % /* update pulse status */
    % *t = *(utsLine + *utsi);
    t = utsLine(utsi+1);
    VVar.t = t;
    t
    % *dt = *(utsLine + (int)fmin(*utsi+1, MaxutsStep-1))-*(utsLine + *utsi);
    dt = utsLine(min(utsi+2, MaxutsStep-1))-utsLine(utsi+1);
    VVar.dt = dt;
    utsi = min(utsi+2, MaxutsStep-1);
    VVar.utsi = utsi;

    if (flag(1) >= 1) % /* update rfAmp, rfPhase, rfFreq, rfCoil for multiple rf lines*/
        for j=0:(flag(1)-1)
            rfCoil = rfCoilLine(rfi+1);
            VVar.rfCoil = rfCoil;
            TxCoili = int32(rfCoil);
            s = rfi + 1;
            while (s < MaxrfStep)
                if (rfCoil == rfCoilLine(s+1))
                    if (abs(rfAmpLine(rfi+1)) <= abs(rfAmpLine(s+1)))
                        rfAmp(TxCoili)= rfAmpLine(rfi+1);
                    else
                        rfAmp(TxCoili)= rfAmpLine(s+1);
                    end
                    if (abs(rfPhaseLine(rfi+1)) <= abs(rfPhaseLine(s+1)))
                        rfPhase(TxCoili)= rfPhaseLine(rfi+1);
                    else
                        rfPhase(TxCoili)= rfPhaseLine(s+1);
                    end
                    if (abs(rfFreqLine(rfi+1)) <= abs(rfFreqLine(s+1)))
                        rfFreq(TxCoili)= rfFreqLine(rfi+1);
                    else
                        rfFreq(TxCoili)= rfFreqLine(s+1);
                    end
                    VVar.rfAmp = rfAmp;
                    VVar.rfFreq = rfFreq;
                    VVar.rfPhase = rfPhase;
                    break;
                end
                s = s + 1;
            end
            
            rfi = rfi + 1;
            VVar.rfi = rfi;
        end
    end
    
    if (flag(2)==1)% /*update GzAmp */
        if (abs(GzAmpLine(Gzi+1)) <= abs(GzAmpLine(int32(min(Gzi+2, MaxGzStep)))))
            GzAmp = GzAmpLine(Gzi+1);
            VVar.GzAmp = GzAmp;
        else
            GzAmp = GzAmpLine(Gzi+2);
            VVar.GzAmp = GzAmp;
        end
        Gzi = Gzi + 1;
        VVar.Gzi = Gzi;
    end
    
    if (flag(3)==1 )% /*update GyAmp */
        if (abs(GyAmpLine(Gyi+1)) <= abs(GyAmpLine(int32(min(Gyi+2, MaxGyStep)))))
            GyAmp = GyAmpLine(Gyi+1);
            VVar.GyAmp = GyAmp;
        else
            GyAmp = GyAmpLine(Gyi+2);
            VVar.GyAmp = GyAmp;
        end
        Gyi = Gyi + 1;
        VVar.Gyi = Gyi;
    end
    
    if (flag(4)==1 )% /*update GxAmp */
        if (abs(GxAmpLine(Gxi+1)) <= abs(GxAmpLine(int32(min(Gxi+2, MaxGxStep)))))
            GxAmp = GxAmpLine(Gxi+1);
            VVar.GxAmp = GxAmp;
        else
            GxAmp = GxAmpLine(Gxi+2);
            VVar.GxAmp = GxAmp;
        end
        Gxi = Gxi + 1;
        VVar.Gyi = Gyi;
    end
    
    ADC = 0;   %/* prevent ADC overflow*/
    VVar.ADC = ADC;
    if (flag(5)==1)% /* update ADC */
        ADC = ADCLine(ADCi+1);
        VVar.ADC = ADC;
        ADCi = ADCi + 1;
        VVar.ADCi = ADCi;
    end
    
    if (flag(6)==1)% /* update Ext */
        Ext = ExtLine(Exti+1);
        VVar.Ext=Ext;
        %/* execute extended process */
        if (Ext ~= 0)
            %ExtCall = mexEvalString("DoExtPlugin");
            %ExtCall = DoExtPlugin;
            disp("DoExtPlugin: ERROR NOT IMPLEMENTED");
            ExtCall = false;
            if (ExtCall)
                disp("Extended process encounters ERROR!");
                return;
            end
            
            % /* update pointers, avoid pointer change between Matlab and Mex call */
            % MzBase          = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "Mz"));
            % MyBase          = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "My"));
            % MxBase          = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "Mx"));
            % RhoBase         = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "Rho"));
            % T1Base          = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "T1"));
            % T2Base          = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VObj"), 0, "T2"));
            % dWRndBase       = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VMag"), 0, "dWRnd"));
            % dB0Base         = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VMag"), 0, "dB0"));
            % GzgridBase      = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VMag"), 0, "Gzgrid"));
            % GygridBase      = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VMag"), 0, "Gygrid"));
            % GxgridBase      = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VMag"), 0, "Gxgrid"));
            MzBase          = VObj.Mz;
            MyBase          = VObj.My;
            MxBase          = VObj.Mx;
            RhoBase         = VObj.Rho;
            T1Base          = VObj.T1;
            T2Base          = VObj.T2;
            dB0Base         = VMag.dB0;
            dWRndBase       = VMag.dWRnd;
            GzgridBase      = VMag.Gzgrid;
            GygridBase      = VMag.Gygrid;
            GxgridBase      = VMag.Gxgrid;
            
            % TxCoilmgBase    = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VCoi"), 0, "TxCoilmg"));
            % TxCoilpeBase    = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VCoi"), 0, "TxCoilpe"));
            % RxCoilx         = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VCoi"), 0, "RxCoilx"));
            % RxCoily         = (float*) mxGetData(mxGetField(mexGetVariablePtr("global", "VCoi"), 0, "RxCoily"));
            TxCoilmgBase    = VCoi.TxCoilmg;
            TxCoilpeBase    = VCoi.TxCoilpe;
            RxCoilx         = VCoi.RxCoilx;
            RxCoily         = VCoi.RxCoily;
            
            % t               = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "t"));
            % dt              = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "dt"));
            % rfAmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "rfAmp"));
            % rfPhase         = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "rfPhase"));
            % rfFreq          = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "rfFreq"));
            % rfCoil          = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "rfCoil"));
            % rfRef           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "rfRef"));
            % GzAmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "GzAmp"));
            % GyAmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "GyAmp"));
            % GxAmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "GxAmp"));
            %                 ADC             = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "ADC"));
            %                 Ext             = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Ext"));
            %                 KzTmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Kz"));
            %                 KyTmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Ky"));
            %                 KxTmp           = (double*) mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Kx"));
            %                 utsi            = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "utsi"));
            %                 rfi             = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "rfi"));
            %                 Gzi             = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Gzi"));
            %                 Gyi             = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Gyi"));
            %                 Gxi             = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Gxi"));
            %                 ADCi            = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "ADCi"));
            %                 Exti            = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "Exti"));
            %                 TRCount         = (int*)    mxGetData(mxGetField(mexGetVariablePtr("global", "VVar"), 0, "TRCount"));
            t               = VVar.t;
            dt              = VVar.dt;
            rfAmp           = VVar.rfAmp;
            rfPhase         = VVar.rfPhase;
            rfFreq          = VVar.rfFreq;
            rfCoil          = VVar.rfCoil;
            rfRef           = VVar.rfRef;
            GzAmp           = VVar.GzAmp;
            GyAmp           = VVar.GyAmp;
            GxAmp           = VVar.GxAmp;
            ADC             = VVar.ADC;
            Ext             = VVar.Ext;
            KzTmp           = VVar.Kz(1);
            KyTmp           = VVar.Ky(1);
            KxTmp           = VVar.Kx(1);
            utsi            = VVar.utsi;
            rfi             = VVar.rfi;
            Gzi             = VVar.Gzi;
            Gyi             = VVar.Gyi;
            Gxi             = VVar.Gxi;
            ADCi            = VVar.ADCi;
            Exti            = VVar.Exti;
            TRCount         = VVar.TRCount;

        end
        Exti = Exti + 1;
        VVar.Exti = Exti;
    end
    
    if (flag(1)+flag(2)+flag(3)+flag(4)+flag(5)+flag(6) == 0)% /* reset VVar */
        rfAmp(1:TxCoilNum)=0;
        rfPhase(1:TxCoilNum)=0;
        rfFreq(1:TxCoilNum)=0;
        %ippsZero_64f(rfAmp, *TxCoilNum);
        %ippsZero_64f(rfPhase, *TxCoilNum);
        %ippsZero_64f(rfFreq, *TxCoilNum);
        GzAmp = 0;
        GyAmp = 0;
        GxAmp = 0;
        ADC = 0;
        Ext = 0;
        %VVar.GzAmp = GzAmp;
        %VVar.GyAmp = GyAmp;
        %VVar.GxAmp = GxAmp;
        VVar.ADC = ADC;
        VVar.Ext = Ext;
    end
    
    %/* execute spin precessing */
    if (dt == 0)% /* end of time point */
        continue;
    elseif (dt < 0)% /* uncontinuous time point process */
        TRCount = TRCount + 1;
        VVar.TRCount = TRCount;
        fprintf("TR Counts: %d of %d\n", TRCount, TRNum);
        
        switch (RunMode)
            case 1 %/* spin rotation simulation & rf simulation */
                disp("SpinRotation: ERROR NOT IMPLEMENTED");
                %MzsBase = MzBase;
                %MysBase = MyBase;
                %MxsBase = MxBase;
                %ippsCopy_32f(MzBase, MzsBase+(*utsi)*(*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum, (*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum);
                %ippsCopy_32f(MyBase, MysBase+(*utsi)*(*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum, (*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum);
                %ippsCopy_32f(MxBase, MxsBase+(*utsi)*(*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum, (*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum);
                Muts(utsi+1) = Muts(utsi);
                break;
        end
        continue;
    end
    
    for Typei = 0:(TypeNum-1)
        for Spini = 0:(SpinNum-1)
            %/* signal acquisition */
            if (ADC == 1)
                % SpinMxNum = num pixels (x,y)
                % SpinMxSliceNum = num slices (3)
                % SpinNum = 1?, Typei = 0?
                %Mx = MxBase + (Typei*(*SpinNum)*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum);
                %My = MyBase + (Typei*(*SpinNum)*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum);
                offset = (Typei*SpinNum*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum+1);
                Mx = MxBase(offset:(offset+SpinMxSliceNum*SpinMxNum-1));
                My = MyBase(offset:(offset+SpinMxSliceNum*SpinMxNum-1));
                for RxCoili = 0:(RxCoilNum-1)%  /* signal acquisition per Rx coil */
                    %/* RxCoil sensitivity */
                    if (RxCoilDefault == 0)
                        disp("RxCoilDefault==0: ERROR NOT IMPLEMENTED");
                        %ippsMul_32f(Mx, RxCoilx+RxCoili*SpinMxSliceNum*SpinMxNum, buffer1, SpinMxSliceNum*SpinMxNum);
                        %ippsMul_32f(My, RxCoily+RxCoili*SpinMxSliceNum*SpinMxNum, buffer2, SpinMxSliceNum*SpinMxNum);
                        %ippsAdd_32f(buffer1, buffer2, buffer3, SpinMxSliceNum*SpinMxNum);
                        
                        %ippsMul_32f(Mx, RxCoily+RxCoili*SpinMxSliceNum*SpinMxNum, buffer1, SpinMxSliceNum*SpinMxNum);
                        %ippsMulC_32f_I(-1, buffer1, SpinMxSliceNum*SpinMxNum);
                        %ippsMul_32f(My, RxCoilx+RxCoili*SpinMxSliceNum*SpinMxNum, buffer2, SpinMxSliceNum*SpinMxNum);
                        %ippsAdd_32f(buffer1, buffer2, buffer4, SpinMxSliceNum*SpinMxNum);
                        
                        %ippsCopy_32f(buffer3, buffer1, SpinMxSliceNum*SpinMxNum);
                        %ippsCopy_32f(buffer4, buffer2, SpinMxSliceNum*SpinMxNum);
                    else
                        %ippsCopy_32f(Mx, buffer1, SpinMxSliceNum*SpinMxNum);
                        %ippsCopy_32f(My, buffer2, SpinMxSliceNum*SpinMxNum);
                        %ippsCopy_32f(Mx, buffer3, SpinMxSliceNum*SpinMxNum);
                        %ippsCopy_32f(My, buffer4, SpinMxSliceNum*SpinMxNum);
                        buffer1 = Mx(1:SpinMxSliceNum*SpinMxNum);
                        buffer2 = My(1:SpinMxSliceNum*SpinMxNum);
                        buffer3 = Mx(1:SpinMxSliceNum*SpinMxNum);
                        buffer4 = My(1:SpinMxSliceNum*SpinMxNum);
                    end
                    
                    %/* rfRef for demodulating rf Phase */
                    if (rfRef ~= 0)
                        disp("rfRef ~= 0: ERROR NOT IMPLEMENTED");
                        %ippsMulC_32f_I((float)cos(-(*rfRef)), buffer1, SpinMxSliceNum*SpinMxNum);
                        %ippsMulC_32f_I((float)-sin(-(*rfRef)), buffer2, SpinMxSliceNum*SpinMxNum);
                        %ippsAdd_32f_I(buffer2, buffer1, SpinMxSliceNum*SpinMxNum);
                        
                        %ippsMulC_32f_I((float)sin(-(*rfRef)), buffer3, SpinMxSliceNum*SpinMxNum);
                        %ippsMulC_32f_I((float)cos(-(*rfRef)), buffer4, SpinMxSliceNum*SpinMxNum);
                        %ippsAdd_32f_I(buffer4, buffer3, SpinMxSliceNum*SpinMxNum);
                    else
                        %ippsCopy_32f(buffer4, buffer3, SpinMxSliceNum*SpinMxNum);
                        buffer3 = buffer4;
                    end
                    
                    %/* signal summation */
                    %ippsSum_32f(buffer1, SpinMxSliceNum*SpinMxNum, buffer, ippAlgHintFast);
                    buffer = sum(buffer1);
                    %Sx[Typei*(*RxCoilNum)*(*SignalNum)+RxCoili*(*SignalNum)+Signali] += (double)*buffer;
                    Sx(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1) = ...
                        Sx(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1) + buffer;
                    VSig.Sx(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1) = ...
                        Sx(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1);
                    %ippsSum_32f(buffer3, SpinMxSliceNum*SpinMxNum, buffer, ippAlgHintFast);
                    buffer = sum(buffer3);
                    %Sy[Typei*(*RxCoilNum)*(*SignalNum)+RxCoili*(*SignalNum)+Signali] += (double)*buffer;
                    Sy(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1) = ...
                        Sy(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1) + buffer;
                    VSig.Sy(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1) = ...
                        Sy(Typei*RxCoilNum*SignalNum+RxCoili*SignalNum+Signali+1);
                end
            end
            
            %/* set openMP for core inner loop */
            
            %#pragma omp parallel num_threads(ThreadNum){
            
            %#pragma omp for private(Slicei, Mz, My, Mx, dWRnd, Rho, T1, T2, dB0, Gzgrid, Gygrid, Gxgrid, TxCoilmg, TxCoilpe)
            %/* need private clause to keep variables isolated, VERY IMPORTANT!!! Otherwise cause cross influence by threads */
            
            for Slicei=0:(SpinMxSliceNum-1)
                
                %/* don't put Matlab function here, may cause Matlab crash! */
                
                %/* set pointer offset for input */
                %Mz			= MzBase + (Typei*(*SpinNum)*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                %My			= MyBase + (Typei*(*SpinNum)*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                %Mx			= MxBase + (Typei*(*SpinNum)*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                %dWRnd		= dWRndBase + (Typei*(*SpinNum)*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                offset = (Typei*SpinNum*SpinMxSliceNum*SpinMxNum + Spini*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum+1);
                Mz			= MzBase(offset:(offset+SpinMxNum-1));
                My			= MyBase(offset:(offset+SpinMxNum-1));
                Mx			= MxBase(offset:(offset+SpinMxNum-1));
                dWRnd		= dWRndBase(offset:(offset+SpinMxNum-1));
                
                %Rho			= RhoBase+(Typei*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                %T1			= T1Base+ (Typei*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                %T2			= T2Base+ (Typei*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum);
                offset = (Typei*SpinMxSliceNum*SpinMxNum + Slicei*SpinMxNum+1);
                Rho			= RhoBase(offset:(offset+SpinMxNum-1));
                T1			= T1Base(offset:(offset+SpinMxNum-1));
                T2			= T2Base(offset:(offset+SpinMxNum-1));
                
                %dB0			= dB0Base + Slicei*SpinMxNum;
                %Gzgrid		= GzgridBase + Slicei*SpinMxNum;
                %Gygrid		= GygridBase + Slicei*SpinMxNum;
                %Gxgrid		= GxgridBase + Slicei*SpinMxNum;
                %TxCoilmg	= TxCoilmgBase + Slicei*SpinMxNum;
                %TxCoilpe	= TxCoilpeBase + Slicei*SpinMxNum;
                offset = Slicei*SpinMxNum+1;
                dB0			= dB0Base(offset:(offset+SpinMxNum-1));
                Gzgrid		= GzgridBase(offset:(offset+SpinMxNum-1));
                Gygrid		= GygridBase(offset:(offset+SpinMxNum-1));
                Gxgrid		= GxgridBase(offset:(offset+SpinMxNum-1));
                TxCoilmg	= TxCoilmgBase(offset:(offset+SpinMxNum-1));
                TxCoilpe	= TxCoilpeBase(offset:(offset+SpinMxNum-1));
                
                % /* call spin discrete precessing */
                try 
                %disp("HERE1");
                [Mx, My, Mz] = BlochKernelNormalCPU(Gyro, CS, SpinNum, Rho, T1, T2, Mz, My, Mx, ...
                    dB0, dWRnd, Gzgrid, Gygrid, Gxgrid, TxCoilmg, TxCoilpe, ...
                    dt, rfAmp, rfPhase, rfFreq, GzAmp, GyAmp, GxAmp, ...
                    Typei, SpinMxNum, SpinMxSliceNum, TxCoilNum);
                %disp("HERE2");
                catch me
                    disp(me.message)
                MxBase(offset:(offset+SpinMxNum-1)) = Mx;
                MyBase(offset:(offset+SpinMxNum-1)) = My;
                MzBase(offset:(offset+SpinMxNum-1)) = Mz;
                if (Slicei==0)
                    figure(2), imshow(abs(reshape(Mx,[100,80])),[]);
                    drawnow;
                end
            end
        end
    end
    
    if (ADC == 1)
        Kz(Signali+1) = Kz(Signali+1) + KzTmp;
        Ky(Signali+1) = Ky(Signali+1) + KyTmp;
        Kx(Signali+1) = Kx(Signali+1) + KxTmp;
        VSig.Kz(Signali+1) = Kz(Signali+1);
        VSig.Ky(Signali+1) = Ky(Signali+1);
        VSig.Kx(Signali+1) = Kx(Signali+1);
        Signali = Signali + 1;
    end
    
    %/* update Kz, Ky & Kx */
    KzTmp = KzTmp + (GzAmp)*(dt)*(Gyro/(2*pi));
    KyTmp = KyTmp + (GyAmp)*(dt)*(Gyro/(2*pi));
    KxTmp = KxTmp + (GxAmp)*(dt)*(Gyro/(2*pi));
    VVar.KzTmp = KzTmp;
    VVar.KyTmp = KyTmp;
    VVar.KxTmp = KxTmp;
    
    switch (RunMode)
        case 1 %/* spin rotation simulation & rf simulation */
            %ippsCopy_32f(MzBase, MzsBase+(*utsi)*(*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum, (*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum);
            %ippsCopy_32f(MyBase, MysBase+(*utsi)*(*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum, (*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum);
            %ippsCopy_32f(MxBase, MxsBase+(*utsi)*(*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum, (*TypeNum)*(*SpinNum)*SpinMxSliceNum*SpinMxNum);
            Muts(utsi+1) = Muts(utsi) + dt;
            break;            
    end
end
%VSig.Mx = Mx;
%VSig.My = My;
%VSig.Mz = Mz;

end
