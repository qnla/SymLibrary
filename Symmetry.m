 classdef Symmetry < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess = 'public', SetAccess = 'private')
        Fusion2 = [];
        Dim = []; %Dim is (number, dimension size, label)
        Name = [];
        TrivialIrrep = [];
        InverseIrrep = [];
        
        F = [];
        FSq = [];
        %FLabels =[];
        %FCumSum = [];
        
        R = [];
        RCW = [];
        RCCW = [];
        B = [];
        Frobenius_Schur_Indicator = [];
        Round = [];
        
        MaxCharges = 0;
        MaxMultiplicities = 0;
    end
    
    properties(GetAccess = 'public', SetAccess = 'private')
        FlagSym = true;
        FlagBraiding = false;
        FlagNonAbelian = false;
        FlagMultiplicities = false;
        FlagRibbon = false;
        FlagSelfDual = false;
        FlagSphericalSymmetry = false;
        
        FlagSymCon_CheckIfEmpty = false;
        FlagSymCon_MaxTensorSize = 2^(31-4);%2GB but in blocks of 16Bytes
        
        
        FlagSaved = false([0,1]);
        SavedInfoDetails = cell([0,23]);
        SavedCheckDetails = cell([0,6]);
        
        FlagSavedMid = false([0,1]);
        SavedCheckMidDetails = cell([0,3]);
        SavedInfoMidDetails = cell([0,9]);
        
        FlagSavedFine = false([0,1]);
        SavedCheckFineDetails = cell([0,1]);
        SavedInfoFineDetails = cell([0,1]);
        
        FlagSavedPlanar = false([0,1]);
        SavedCheckPlanarDetails = cell([0,0]);
        SavedInfoPlanarDetails = cell([0,8]);
        
        FlagSavedMidPlanar = false([0,1]);
        %SavedCheckMidPlanarDetails = cell([0,0]);
        %SavedInfoMidPlanarDetails = cell([0,5]);
    end
    
    
    %general Symmetry functions
    methods
        
        function self = FixUp(self)
            self.FlagSymCon_CheckIfEmpty = false;
        end
        
        function BoolVal = IsBraiding(self)
            BoolVal = self.FlagBraiding;
        end

        function BoolVal = IsNonAbelian(self)
            BoolVal = self.FlagNonAbelian;
        end

        function BoolVal = IsSelfDual(self)
            BoolVal = self.FlagSelfDual;
        end

        function BoolVal = IsMultiplicities(self)
            BoolVal = self.FlagMultiplicities;
        end
        
        function BoolVal = IsSpherical(self)
            BoolVal = self.FlagSphericalSymmetry;
        end
        
        function BoolVal = IsRibbon(self)
            BoolVal = self.FlagRibbon;
        end
        
        function BoolVal = IsSym(self)
            BoolVal = self.FlagSym;
        end
        
        function Out = getDim(self)
            Out = self.Dim();
        end
        
        function Out = getDimExternal(self)
            Out = self.Dim();
        end
        
        function Out = getFusion2(self)
            Out = self.Fusion2();
        end
        
        function Out = getInverseIrrep(self)
            Out = self.InverseIrrep();
        end
        
        function Out = getTrivialIrrep(self)
            Out = self.TrivialIrrep();
        end
        
        function [Internal, MultList, Multiplicities] = FuseChargeList(self,A,B)
            Fuse = reshape(self.Fusion2, [size(self.Fusion2,1),size(self.Fusion2,2)*size(self.Fusion2,3)]);
            Results = Fuse(:,B+(A-1)*size(self.Fusion2,3))';
            Temp = cumsum(sum(Results',1));
            
            %this is a list which allows us to extend our External charges
            %to account for multiple charges.
            MultList = sum(repmat((1:Temp(end))',[1,size(Temp,2)])>repmat(Temp,[Temp(end),1]),2)+1;
            
            Results = Results';
            TempSing = cumsum(Results(:)');
            if TempSing(end)~=Temp(end)
                error('Affirmation Failed: There are different numbers of outputs if I do this singly or as a group');
            end
            
            %this is a list of what the internal charges will take, however
            %it only occurs once, so if there are multiplicities we need to
            %multiply the number of entries.
            TempValues = mod(find(Results(:))-1,size(Results,1))+1;
            
            %This should be the number of multiplicities for each fusion
            %which we need to work out for extending TempValues to Internal:
            TempSingMod = TempSing(Results(:)~=0);
            Multiplicities = sum(repmat((1:TempSing(end))',[1,size(TempSingMod,2)])>repmat(TempSingMod,[TempSing(end),1]),2)+1;
            Results = Results';
            
            %values for the Internal charges that come out (multiplying by
            %the number of multiplicities:
            Internal = TempValues(Multiplicities);
            
            %The actual Multiplicities Counter
            Multiplicities= (1:TempSing(end))'-sum(repmat(TempSing(1:(end-1)),[TempSing(end),1]).*...
                (((repmat((1:TempSing(end))',[1,size(TempSing,2)-1])>repmat(TempSing(1:(end-1)),[TempSing(end),1]))    )&...
                (repmat((1:TempSing(end))',[1,size(TempSing,2)-1])<=repmat(TempSing(2:(end)),[TempSing(end),1]))),2);
        end
        
        
        function [OutF, OutFCumSum, OutFLabels] = getFMoves(self)
            OutF = self.F;
            OutFCumSum = self.FCumSum;
            OutFLabels = self.FLabels;
        end
        
        function OutR = getRMoves(self)
            OutR = self.R;
        end
        
        function OutName = getSymName(self)
            OutName = self.Name;
        end
        
        function [OutCharges, OutNames] = FuseCharges(self, a, b, aNames, bNames)
            N = self.Fusion2;
            Dim = self.getDim;
            
            if length(size(a)) ~= 2
                error('size of first input is wrong')
            elseif ~any(size(a) == 1)
                error('size of first input is wrong')
            end
            
            if length(size(b)) ~= 2
                error('size of first input is wrong')
            elseif ~any(size(b) == 1)
                error('size of first input is wrong')
            end
            
            a = a(:);
            b = b(:);
            
            if nargin < 5
                KeepB = true(size(b));
            else
                if length(size(bNames)) ~= 2
                    error('names for first input is wrong')
                elseif size(bNames,2) ~= size(Dim,2)-2
                    error('names for first input is wrong, needs to be a column of lables (which are lists of numbers)')
                end
                
                KeepB = any(all(repmat(reshape(Dim(:,3:end), [size(Dim,1),1,size(Dim,2)-2]),[1,size(bNames,1)])...
                              == repmat(reshape(bNames, [1, size(bNames,1), size(bNames,2)]),[size(Dim,1),1,1]),3),2);
                
            end
            
            if nargin < 4
                KeepA = true(size(a));
            else
                if length(size(aNames)) ~= 2
                    error('names for first input is wrong')
                elseif size(aNames,2) ~= size(Dim,2)-2
                    error('names for first input is wrong, needs to be a column of lables (which are lists of numbers)')
                end
                
                KeepA = any(all(repmat(reshape(Dim(:,3:end), [size(Dim,1),1,size(Dim,2)-2]),[1,size(aNames,1)])...
                              == repmat(reshape(aNames, [1, size(aNames,1), size(aNames,2)]),[size(Dim,1),1,1]),3),2);
                
            end
            
            if all(~KeepA) | all(~KeepB)
                OutCharges = [];
                OutNames = [];
                return;
            end
            
            N = N(:,KeepA,KeepB);
            
            SizeN = size(N);
            if size(N)<3
                SizeN = [SizeN,ones([1,3-size(SizeN,2)])];
            end 
            N = reshape(N, [SizeN(1)*SizeN(2), SizeN(3)])*b;
            N = reshape(N, [SizeN(1),SizeN(2)])*a;
            
            KeepOut = N~=0;
            
            if nargout>1
                OutCharges = N(KeepOut);
                OutNames = Dim(KeepOut,3:end);
            else
                OutCharges = N;
            end
            
        end
        
        function OutputConvertNumbers = GenerateNumbers(self, WordDescriptions, ChargesExternal, ChargesInternal, Multiplicities)
            
            if ~iscell(WordDescriptions)
                error('Error: This can''t be the correct Descriptions for the processes')
            end
            
            
            if isempty(WordDescriptions)
                EntriesS = ones([size(ChargesExternal,1),1]);
                ListIn = 1:size(ChargesExternal,1);
                ListOut = 1:size(ChargesExternal,1);
                
                LabelsInAll = [ChargesExternal, ChargesInternal, Multiplicities];
                LabelsOut = [ChargesExternal, ChargesInternal, Multiplicities];
                
                OutputConvertNumbers = {LabelsInAll,LabelsOut,sparse(ListIn,ListOut,EntriesS,max(ListIn),max(ListOut))};
                return;
            end
            
            if ~iscell(WordDescriptions{1})
                FlagOutputSingle = true;
                WordDescriptions = {WordDescriptions};
                ChargesExternal = {ChargesExternal};
                ChargesInternal = {ChargesInternal};
                Multiplicities = {Multiplicities};
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Still need to check that this makes sense here.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            FlagDoNothing = false;
            FlagIgnoreBraiding = false;
            FlagIgnoreFmoves = false;
            if ~self.FlagBraiding
                if ~self.FlagNonAbelian
                    FlagDoNothing = true;
                else %then only non-abelian
                    FlagIgnoreBraiding = true;
                end
            else
                if ~self.FlagNonAbelian
                    FlagIgnoreFmoves = true;
                else %then only non-abelian
                    %we do everything here then.
                end
            end
            
            
            ChargeMax = self.MaxCharges;
            MultMax = self.MaxMultiplicities;
            for kk = 1:numel(WordDescriptions)
                
                if isempty(WordDescriptions{kk})
                    EntriesS = ones([size(ChargesExternal{kk},1),1]);
                    ListIn = 1:size(ChargesExternal{kk},1);
                    ListOut = 1:size(ChargesExternal,1);
                    
                    LabelsInAll = [ChargesExternal{kk}, ChargesInternal{kk}, Multiplicities{kk}];
                    LabelsOut = [ChargesExternal{kk}, ChargesInternal{kk}, Multiplicities{kk}];
                    
                    OutputConvertNumbers{kk} = {LabelsInAll,LabelsOut,sparse(ListIn,ListOut,EntriesS,max(ListIn),max(ListOut))};
                    continue;
                end
                InitialExtNumber = size(ChargesExternal{kk},2);
                InitialIntNumber = size(ChargesInternal{kk},2);
                InitialMultNumber = size(Multiplicities{kk},2);
                FinalExtNumber = size(ChargesExternal{kk},2);
                FinalIntNumber = size(ChargesInternal{kk},2);
                FinalMultNumber = size(Multiplicities{kk},2);
                
                EntriesS = ones([size(ChargesExternal{kk},1),1]);
                LabelsInAll = [ChargesExternal{kk}, ChargesInternal{kk}, Multiplicities{kk}];
                Labels = [LabelsInAll,LabelsInAll];
                %the first 3 are our initial charges, the last 3 are our
                %final charges, these will change with iterations.
                
                for ll = 1:size(WordDescriptions{kk},1)
                    if ~isequal(WordDescriptions{kk}{ll,1}, 'Swap')&&~isequal(WordDescriptions{kk}{ll,1}, 'ForceSame')
                        FlagSkip = false;
                        
                        switch WordDescriptions{kk}{ll,1}
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %Fusion Actions
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            case 'F'
                                
                                C = 6;
                                
                                LabelsWorking = self.F{1};
                                Entries = self.F{2};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3,4,6];
                                OutputCharges = [1,2,3,4,5];
                                if self.IsMultiplicities %check these are correct
                                    InputMult = [9,10];
                                    OutputMult = [7,8];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'FInv'
                                
                                C = 6;
                                LabelsWorking = self.F{1};
                                Entries = self.F{3};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3,4,5];
                                OutputCharges = [1,2,3,4,6];
                                if self.IsMultiplicities
                                    OutputMult = [9,10];
                                    InputMult = [7,8];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'FSq'
                                
                                C = 6;
                                LabelsWorking = self.FSq{1};
                                Entries = self.FSq{2};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3,4,6];
                                OutputCharges = [1,2,3,4,5];
                                if self.IsMultiplicities
                                    InputMult = [9,10];
                                    OutputMult = [7,8];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'FSqInv'
                                
                                C = 6;
                                LabelsWorking = self.FSq{1};
                                Entries = self.FSq{3};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3,4,5];
                                OutputCharges = [1,2,3,4,6];
                                if self.IsMultiplicities
                                    OutputMult = [9,10];
                                    InputMult = [7,8];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                                
                                
                                
                                
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %Braiding Actions
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                            case 'R'
                                
                                C = 3;
                                LabelsWorking = self.R{1};
                                Entries = self.R{2};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3];
                                OutputCharges = [2,1,3];
                                if self.IsMultiplicities
                                    InputMult = [5];
                                    OutputMult = [4];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                [~,Keep,~] = unique(Details(1,InputCharges));
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                if sum(InputCharges==1|InputCharges==2)==2
                                    Loc1 = InputCharges == 1;
                                    Loc2 = InputCharges == 2;
                                    OutputCharges = InputCharges;
                                    OutputCharges(Loc1)=2;
                                    OutputCharges(Loc2)=1;
                                else
                                    OutputCharges = InputCharges;
                                end
                                
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                OutputCharges = Details(1,InputCharges);
                                InputCharges = Details(1,InputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'RInv'
                                
                                C = 3;
                                LabelsWorking = self.R{1};
                                Entries = self.R{3};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3];
                                OutputCharges = [2,1,3];
                                if self.IsMultiplicities
                                    OutputMult = [5];
                                    InOputMult = [4];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                [~,Keep,~] = unique(Details(1,InputCharges));
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                if sum(InputCharges==1|InputCharges==2)==2
                                    Loc1 = InputCharges == 1;
                                    Loc2 = InputCharges == 2;
                                    OutputCharges = InputCharges;
                                    OutputCharges(Loc1)=2;
                                    OutputCharges(Loc2)=1;
                                else
                                    OutputCharges = InputCharges;
                                end
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                OutputCharges = Details(1,InputCharges);
                                InputCharges = Details(1,InputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'B'
                                
                                C = 6;
                                LabelsWorking = self.B{1};
                                Entries = self.B{2};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3,4,6];
                                OutputCharges = [1,3,2,4,5];
                                if self.IsMultiplicities
                                    InputMult = [9,10];
                                    OutputMult = [7,8];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                OutputCharges = Details(1,InputCharges.*(InputCharges<5)+5*(InputCharges == 6));%OutputCharges);
                                InputCharges = Details(1,InputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'BInv'
                                
                                C = 6;
                                LabelsWorking = self.B{1};
                                Entries = self.B{3};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3,4,6];
                                OutputCharges = [1,3,2,4,5];
                                if self.IsMultiplicities
                                    OutputMult = [9,10];
                                    InputMult = [7,8];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                OutputCharges = Details(1,InputCharges.*(InputCharges<5)+5*(InputCharges == 6));%OutputCharges);
                                InputCharges = Details(1,InputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                                
                                
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %Bending Actions
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                            case 'RCW'
                                
                                C = 3;
                                LabelsWorking = self.RCW{1};
                                Entries = self.RCW{2};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3];
                                OutputCharges = [1,2,3];
                                if self.IsMultiplicities
                                    InputMult = [5];
                                    OutputMult = [4];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                [~,Keep,~] = unique(Details(1,InputCharges));
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                OutputCharges = OutputCharges(Keep);
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'RCWInv'
                                
                                C = 3;
                                LabelsWorking = self.RCW{1};
                                Entries = self.RCW{3};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3];
                                OutputCharges = [1,2,3];
                                if self.IsMultiplicities
                                    OutputMult = [5];
                                    InputMult = [4];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                [~,Keep,~] = unique(Details(1,InputCharges));
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                OutputCharges = OutputCharges(Keep);
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    %OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    %OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'RCCW'
                                
                                C = 3;
                                LabelsWorking = self.RCCW{1};
                                Entries = self.RCCW{2};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3];
                                OutputCharges = [1,2,3];
                                if self.IsMultiplicities
                                    InputMult = [5];
                                    OutputMult = [4];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                [~,Keep,~] = unique(Details(1,InputCharges));
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                OutputCharges = OutputCharges(Keep);
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    %InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    %InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end
                                
                            case 'RCCWInv'
                                
                                C = 3;
                                LabelsWorking = self.RCCW{1};
                                Entries = self.RCCW{3};
                                Details = WordDescriptions{kk}{ll,2};
                                if ~self.IsMultiplicities
                                    LabelsWorking = LabelsWorking(:,1:C); %just in case
                                    Details = Details(:,1:C);
                                end
                                
                                ThingsWeCanIgnore = find(isnan(Details(1,:)));
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    Keep = all(LabelsWorking(:,ThingsWeCanIgnore) == repmat((ThingsWeCanIgnore>C)*1+(ThingsWeCanIgnore<=C)*self.TrivialIrrep,[size(LabelsWorking,1),1]),2);
                                    LabelsWorking = LabelsWorking(Keep,:);
                                    Entries = Entries(Keep,:);
                                end
                                
                                ChargesToInvert = find(Details(2,1:C) == -1);
                                LabelsWorking(:,ChargesToInvert) = self.InverseIrrep(LabelsWorking(:,ChargesToInvert));
                                
                                InputCharges = [1,2,3];
                                OutputCharges = [1,2,3];
                                if self.IsMultiplicities
                                    OutputMult = [5];
                                    InputMult = [4];
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end 
                                
                                if ~isempty(ThingsWeCanIgnore)
                                    InputCharges(any(repmat(InputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputCharges,1),1]),1)) = [];
                                    OutputCharges(any(repmat(OutputCharges,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputCharges,1),1]),1)) = [];
                                    InputMult(any(repmat(InputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(InputMult,1),1]),1)) = [];
                                    OutputMult(any(repmat(OutputMult,[1,numel(ThingsWeCanIgnore)])==repmat(ThingsWeCanIgnore(:),[size(OutputMult,1),1]),1)) = [];
                                    
                                    
                                    if isempty(InputCharges)
                                        InputCharges = zeros([1,0]);
                                    end
                                    if isempty(OutputCharges)
                                        OutputCharges = zeros([1,0]);
                                    end
                                    if isempty(InputMult)
                                        InputMult = zeros([1,0]);
                                    end
                                    if isempty(OutputMult)
                                        OutputMult = zeros([1,0]);
                                    end
                                    
                                end
                                
                                [~,Keep,~] = unique(Details(1,InputCharges));
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                OutputCharges = OutputCharges(Keep);
                                
                                if ~self.IsMultiplicities
                                    InputMax = self.MaxCharges^(numel(InputCharges));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,InputCharges);
                                        %InputValues = sum((LabelsWorking(:,InputCharges)-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    OutputMax = self.MaxCharges^(numel(OutputCharges));
                                    OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,OutputCharges);
                                        %OutputValues = sum((LabelsWorking(:,OutputCharges)-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                else
                                    InputMax = self.MaxCharges^(numel(InputCharges))*self.MaxMultiplicities^(numel(InputMult));
                                    InputSizes = [self.MaxCharges.^(0:(numel(InputCharges)-1)),self.MaxCharges.^(numel(InputCharges))*self.MaxMultiplicities^(0:(numel(InputMult)-1))];
                                    if ~isempty(InputCharges)
                                        LeftB = LabelsWorking(:,[InputCharges,InputMult]);
                                        %InputValues = sum((LabelsWorking(:,[InputCharges,InputMult])-1).*repmat(InputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelsWorking,1),1]);
                                        %InputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                    
                                    
                                    %OutputMax = self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities^(numel(OutputMult));
                                    %OutputSizes = [self.MaxCharges.^(0:(numel(OutputCharges)-1)),self.MaxCharges^(numel(OutputCharges))*self.MaxMultiplicities.^(0:(numel(OutputMult)-1))];
                                    if ~isempty(OutputCharges)
                                        RightB = LabelsWorking(:,[OutputCharges,OutputMult]);
                                        %OutputValues = sum((LabelsWorking(:,[OutputCharges,OutputMult])-1).*repmat(OutputSizes,[size(LabelsWorking,1),1]),2)+1;
                                    else
                                        RightB = ones([size(LabelsWorking,1),1]);
                                        %OutputValues = ones([size(LabelsWorking,1),1]);
                                    end
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                                InputCharges = Details(1,InputCharges);
                                OutputCharges = Details(1,OutputCharges);
                                if self.IsMultiplicities
                                    InputMult = Details(1,InputMult);
                                    OutputMult = Details(1,OutputMult);
                                else
                                    InputMult = zeros([1,0]);
                                    OutputMult = zeros([1,0]);
                                end                           
                                
                                
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %phase Actions
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                            case 'Round'
                                InputMult = zeros([1,0]);
                                OutputMult = zeros([1,0]);
                                
                                Details = WordDescriptions{kk}{ll,2};
                                if isnan(Details(1,:))
                                    FlagSkip = true;
                                else
                                    InputCharges = Details(1,1);
                                    OutputCharges = Details(1,1);
                                    if Details(2,1) == -1
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Round(self.InverseIrrep(:));
                                        %Moves = sparse(diag(self.Round(self.InverseIrrep)));
                                    else
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Round(:);
                                        %Moves = sparse(diag(self.Round));
                                    end
                                end
                                
                            case 'RoundInv'
                                InputMult = zeros([1,0]);
                                OutputMult = zeros([1,0]);
                                
                                Details = WordDescriptions{kk}{ll,2};
                                if isnan(Details(1,:))
                                    FlagSkip = true;
                                else
                                    InputCharges = Details(1,1);
                                    OutputCharges = Details(1,1);
                                    if Details(2,1) == -1
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Round(self.InverseIrrep(:)).'';
                                        %Moves = sparse(diag(self.Round(self.InverseIrrep).''));
                                    else
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Round(:).'';
                                        %Moves = sparse(diag(self.Round.''));
                                    end
                                end
                                
                            case 'pDown'
                                InputMult = zeros([1,0]);
                                OutputMult = zeros([1,0]);
                                
                                Details = WordDescriptions{kk}{ll,2};
                                if isnan(Details(1,:))
                                    FlagSkip = true;
                                else
                                    InputCharges = Details(1,1);
                                    OutputCharges = Details(1,1);
                                    if Details(2,1) == -1
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:);
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator));
                                    else
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:).';
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator.'));
                                    end
                                end
                                
                            case 'pDownInv'
                                InputMult = zeros([1,0]);
                                OutputMult = zeros([1,0]);
                                
                                Details = WordDescriptions{kk}{ll,2};
                                if isnan(Details(1,:))
                                    FlagSkip = true;
                                else
                                    InputCharges = Details(1,1);
                                    OutputCharges = Details(1,1);
                                    if Details(2,1) == -1
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:).';
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator.'));
                                    else
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:);
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator));
                                    end
                                end
                                
                            case 'pUp'
                                InputMult = zeros([1,0]);
                                OutputMult = zeros([1,0]);
                                
                                Details = WordDescriptions{kk}{ll,2};
                                if isnan(Details(1,:))
                                    FlagSkip = true;
                                else
                                    InputCharges = Details(1,1);
                                    OutputCharges = Details(1,1);
                                    if Details(2,1) == -1
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:).';
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator.'));
                                    else
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:);
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator));
                                    end
                                end
                                
                            case 'pUpInv'
                                InputMult = zeros([1,0]);
                                OutputMult = zeros([1,0]);
                                
                                Details = WordDescriptions{kk}{ll,2};
                                if isnan(Details(1,:))
                                    FlagSkip = true;
                                else
                                    InputCharges = Details(1,1);
                                    OutputCharges = Details(1,1);
                                    if Details(2,1) == -1
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:);
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator));
                                    else
                                        LeftB = (1:numel(self.Frobenius_Schur_Indicator))';
                                        RightB = LeftB;
                                        EntriesB = self.Frobenius_Schur_Indicator(:).';
                                        %Moves = sparse(diag(self.Frobenius_Schur_Indicator.'));
                                    end
                                end
                                
                                
                                
                                
                                
                                
                                
                                
                            case 'Dim'
                                
                                NumberChargesCombined = size(WordDescriptions{kk}{ll,2},2)-1;
                                
                                %work out the values of different things:
                                
                                InputChargesTop = [WordDescriptions{kk}{ll,2}(1,2:end),WordDescriptions{kk}{ll,2}(3,2:end)];
                                InputChargesPrime = [WordDescriptions{kk}{ll,2}(2,2:end),WordDescriptions{kk}{ll,2}(4,2:end)];
                                InputChargesMid = WordDescriptions{kk}{ll,2}(1,1);
                                InputMult = WordDescriptions{kk}{ll,2}(5,2:end);
                                InputMultPrime = WordDescriptions{kk}{ll,2}(6,2:end);
                                
                                OutputCharges = WordDescriptions{kk}{ll,2}(3,end);%,WordDescriptions{kk}{ll,2}(4,end)];
                                
                                InputDirections = [WordDescriptions{kk}{ll,2}(1+6,2:end),WordDescriptions{kk}{ll,2}(3+6,2:end)];
                                InputDirectionsPrime = [WordDescriptions{kk}{ll,2}(2+6,2:end),WordDescriptions{kk}{ll,2}(4+6,2:end)];
                                InputDirectionsMid = WordDescriptions{kk}{ll,2}(1+6,1);
                                
                                %InputDirections = [WordDescriptions{kk}{ll,2}(1,2:(end-1)),WordDescriptions{kk}{ll,2}(3,end)];
                                %InputDirectionsPrime = [WordDescriptions{kk}{ll,2}(1,2:(end-1)),WordDescriptions{kk}{ll,2}(4,end)];
                                %InputDirectionsMid = WordDescriptions{kk}{ll,2}(1,1);
                                
                                InputDirectionsAll = InputDirections(1:NumberChargesCombined).*InputDirectionsPrime(1:NumberChargesCombined);
                                
                                
                                %now work out the Labels for the charges
                                
                                InputChargesDim = [InputChargesTop,InputChargesMid,InputChargesPrime];
                                
                                InputDimCharges = InitialExtNumber + InitialIntNumber + InitialMultNumber + abs(InputChargesDim) + (InputChargesDim>0).*FinalExtNumber;
                                
                                InputDimMult = InitialExtNumber + InitialIntNumber + InitialMultNumber + FinalExtNumber + FinalIntNumber + [InputMult,InputMultPrime];
                                
                                LabelCharges = repmat(self.TrivialIrrep, [size(Labels,1),numel(InputDimCharges)]);
                                LabelCharges(:,~isnan(InputDimCharges)) = Labels(:,InputDimCharges(~isnan(InputDimCharges)));
                                
                                LabelMults = ones([size(Labels,1),numel(InputDimMult)]);
                                LabelMults(:,~isnan(InputDimMult)) = Labels(:,InputDimMult(~isnan(InputDimMult)));
                                
                                [~,Keep,~] = unique([LabelCharges,LabelMults],'rows');
                                %[~,Keep,~] = winunique([LabelCharges,LabelMults]);
                                
                                LabelCharges = LabelCharges(Keep,:);
                                LabelMults = LabelMults(Keep,:);
                                
                                
                                
                                %we product the charges of the negative
                                %external charges which are traced over
                                % + the middle charge (MidCharge), then 
                                %divide by the last internal charge
                                
                                Entries = sqrt(prod(reshape(self.Dim(LabelCharges(:,[(1:NumberChargesCombined),numel(InputChargesTop)+1]),2),...
                                    [size(LabelCharges,1),NumberChargesCombined+1]),2)...
                                    ./self.Dim(LabelCharges(:,numel(InputChargesTop)),2));
                                
                                %this will adjust for the direction of
                                %external legs which we are tracing over
                                InputDirectionsPrime(1:NumberChargesCombined) = -InputDirectionsPrime(1:NumberChargesCombined);
                                
                                InputDirectionsAll = [InputDirections,InputDirectionsMid,InputDirectionsPrime];
                                
                                
                                for nn = 1:numel(InputDirectionsAll)
                                    if InputDirectionsAll(nn) == -1 %wrong way (this is why we swapped on Directions Prime earlier
                                        LabelCharges(:,nn) = self.InverseIrrep(LabelCharges(:,nn));
                                    end
                                end
                                
                                
                                Keep2 = all(LabelCharges(:,1:(2*NumberChargesCombined))==LabelCharges(:,2*NumberChargesCombined+1+(1:(2*NumberChargesCombined))),2)...
                                    &all(LabelMults(:,1:NumberChargesCombined)==LabelMults(:,NumberChargesCombined+(1:NumberChargesCombined)),2);
                                
                                LabelCharges = LabelCharges(Keep2,:);
                                LabelMults = LabelMults(Keep2,:);
                                
                                Entries = Entries(Keep2);
                                
                                
                                %do some stuff with trivals:
                                
                                
                                InputCharges = [InputChargesTop, InputChargesMid, InputChargesPrime];
                                InputMult = [InputMult,InputMultPrime];
                                
                                ThingsToForgetCharges = isnan(InputCharges);%(1:(2*NumberChargesCombined+1)));
                                %ThingsToForgetChargesPrime = isnan(InputCharges((2*NumberChargesCombined+2):end));
                                
                                if ThingsToForgetCharges(2*NumberChargesCombined)
                                    FlagTrivialTop = true;
                                else
                                    FlagTrivialTop = false;
                                end
                                
                                if ThingsToForgetCharges(end)
                                    FlagTrivialBottom = true;
                                else
                                    FlagTrivialBottom = false;
                                end
                                
                                %ThingsToForgetCharges(1:(2*NumberChargesCombined)) = ThingsToForgetCharges(1:(2*NumberChargesCombined))|ThingsToForgetChargesPrime(1:(2*NumberChargesCombined));
                                ThingsToForgetMults = isnan(InputMult);%(1:(NumberChargesCombined)))|isnan(InputMult((NumberChargesCombined+1):end));
                                
                                for nn = numel(ThingsToForgetCharges):-1:1
                                    if ThingsToForgetCharges(nn)
                                        InputCharges(nn) = [];%+[0,length(InputCharges)-2*NumberChargesCombined]
                                        Temp = any(LabelCharges(:,nn) ~= self.TrivialIrrep,2);%+[0,size(LabelCharges,2)-2*NumberChargesCombined]
                                        LabelCharges(Temp,:) = [];
                                        LabelCharges(:,nn) = [];%+[0,size(LabelCharges,2)-2*NumberChargesCombined]
                                        LabelMults(Temp,:) = [];
                                        Entries(Temp) = [];
                                    end
                                end
                                
                                %if self.IsMultiplicities
                                    for nn = numel(ThingsToForgetMults):-1:1
                                        if ThingsToForgetMults(nn)
                                            InputMult(nn) = [];%+[0,length(InputMult)-NumberChargesCombined]
                                            Temp = any(LabelMults(:,nn) ~= 1,2);%+[0,size(LabelMults,2)-NumberChargesCombined]
                                            LabelCharges(Temp,:) = [];
                                            LabelMults(Temp,:) = [];
                                            LabelMults(:,nn) = [];%+[0,size(LabelMults,2)-NumberChargesCombined]
                                            Entries(Temp) = [];
                                        end
                                    end
                                %end
                                
                                
                                
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %now do some edge cases
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                
                                
                                if FlagTrivialTop && FlagTrivialBottom
                                    %in this case we know the output will
                                    %be a number and so don't need to worry
                                    %about it:
                                    OutputCharges = zeros([1,0]);
                                    LabelChargesOutput = zeros([size(LabelCharges,1),0]);
                                    OutputMax = 1;
                                    OutputSizes = 1;
                                    
                                    
                                elseif FlagTrivialTop
                                    
                                    %if the bottom is an external leg then
                                    %its easy
                                    if InputChargesPrime(2*NumberChargesCombined)<0
                                        
                                        OutputCharges = InputChargesPrime(2*NumberChargesCombined);
                                        %note that this is still trivial
                                        OutputMax = size(self.Dim,1);
                                        OutputSizes = 1;
                                        LabelChargesOutput = repmat(self.TrivialIrrep,[size(LabelCharges,1),1]);
                                        
                                    elseif all(isnan(WordDescriptions{kk}{ll,2}([2,3],1)))
                                        
                                        %this part is if the external
                                        %charge and internal charges are 
                                        %both NaNs
                                        
                                        OutputCharges = zero([1,0]);
                                        %note that this is still trivial
                                        OutputMax = 1;
                                        OutputSizes = 1;
                                        LabelChargesOutput = repmat(self.TrivialIrrep,[size(LabelCharges,1),1]);
                                        
                                    elseif isnan(WordDescriptions{kk}{ll,2}(3,1))
                                        
                                        %this part is if the external
                                        %charge is an nan as well.
                                        
                                        OutputCharges = WordDescriptions{kk}{ll,2}(2,1);
                                        %note that this is still trivial
                                        OutputMax = size(self.Dim,1);
                                        OutputSizes = 1;
                                        LabelChargesOutput = repmat(self.TrivialIrrep,[size(LabelCharges,1),1]);
                                        
                                    else
                                        %otherwise we need to add in other
                                        %legs:
                                        ExternalExtraLeg = WordDescriptions{kk}{ll,2}(2,1);
                                        InternalExtraLeg = WordDescriptions{kk}{ll,2}(3,1);
                                        ExternalExtraChargeDirection = WordDescriptions{kk}{ll,2}(2+6,1);
                                        InternalExtraChargeDirection = WordDescriptions{kk}{ll,2}(3+6,1);
                                        
                                        %LabelCharges = reshape(LabelCharges,[size(LabelCharges,1),1,size(LabelCharges,2)]);
                                        
                                        LabelChargesOutput = reshape(repmat(reshape(1:size(self.Dim,1),[1,size(self.Dim,1)]),...
                                            [size(LabelCharges,1),1]),[size(self.Dim,1)*size(LabelCharges,1),1]);
                                        
                                        
                                        LabelCharges = [repmat(LabelCharges,[size(self.Dim,1),1]), reshape(repmat(reshape(1:size(self.Dim,1),[1,size(self.Dim,1)]),...
                                            [size(LabelCharges,1),1]),[size(self.Dim,1)*size(LabelCharges,1),1]), reshape(repmat(reshape(1:size(self.Dim,1),[1,size(self.Dim,1)]),...
                                            [size(LabelCharges,1),1]),[size(self.Dim,1)*size(LabelCharges,1),1])];
                                        
                                        
                                        if ExternalExtraChargeDirection ~= +1
                                            LabelCharges(:,end-1) = self.InverseIrrep(LabelCharges(:,end-1));
                                        end
                                        if InternalExtraChargeDirection ~= +1
                                            LabelCharges(:,end) = self.InverseIrrep(LabelCharges(:,end));
                                        end
                                        
                                        
                                        InputCharges = [InputCharges, ExternalExtraLeg,InternalExtraLeg];
                                        OutputCharges = ExternalExtraLeg;
                                        
                                        
                                        Entries = repmat(Entries,[size(self.Dim,1),1]);
                                        
                                        %if self.IsMultiplicities
                                            LabelMults = [repmat(LabelMults,[size(self.Dim,1),1]), ones([size(self.Dim,1)*size(LabelMults,1),1])];
                                            
                                            %will only have one extra mult
                                            %value, as this must be the
                                            %trivial one
                                            
                                            MultiplicitiesExtra = WordDescriptions{kk}{ll,2}(4,1);
                                            InputMult = [InputMult, MultiplicitiesExtra];
                                            
                                        %end
                                        
                                        LabelChargesOutput = LabelCharges(:,end-1);
                                        OutputMax = self.MaxCharges.^(size(LabelChargesOutput,2));
                                        OutputSizes = self.MaxCharges.^(0:(size(LabelChargesOutput,2)-1));
                                        
                                    end
                                    
                                elseif FlagTrivialBottom
                                    %repeat prior case:
                                    
                                    %if the top is an external leg then
                                    %its easy
                                    if InputChargesTop(2*NumberChargesCombined)<0
                                        
                                        LabelChargesOutput = InputChargeTop(2*NumberChargesCombined);
                                        %note that this is still trivial
                                        OutputMax = size(self.Dim,1);
                                        OutputSizes = 1;
                                        LabelChargesOutput = repmat(self.TrivialIrrep,[size(LabelCharges,1),1]);
                                        
                                    elseif all(isnan(WordDescriptions{kk}{ll,2}([2,3],1)))
                                        
                                        %this part is if the external
                                        %charge and internal charges are 
                                        %both NaNs
                                        
                                        OutputCharges = zero([1,0]);
                                        %note that this is still trivial
                                        OutputMax = 1;
                                        OutputSizes = 1;
                                        LabelChargesOutput = repmat(self.TrivialIrrep,[size(LabelCharges,1),1]);
                                        
                                    elseif isnan(WordDescriptions{kk}{ll,2}(3,1))
                                        
                                        %this part is if the external
                                        %charge is an nan as well.
                                        
                                        OutputCharges = WordDescriptions{kk}{ll,2}(2,1);
                                        %note that this is still trivial
                                        OutputMax = size(self.Dim,1);
                                        OutputSizes = 1;
                                        LabelChargesOutput = repmat(self.TrivialIrrep,[size(LabelCharges,1),1]);
                                        
                                    else
                                        %otherwise we need to add in other
                                        %legs:
                                        ExternalExtraLeg = WordDescriptions{kk}{ll,2}(2,1);
                                        InternalExtraLeg = WordDescriptions{kk}{ll,2}(3,1);
                                        ExternalExtraChargeDirection = WordDescriptions{kk}{ll,2}(2+6,1);
                                        InternalExtraChargeDirection = WordDescriptions{kk}{ll,2}(3+6,1);
                                        
                                        %LabelCharges = reshape(LabelCharges,[size(LabelCharges,1),1,size(LabelCharges,2)]);
                                        
                                        LabelChargesOutput = reshape(repmat(reshape(1:size(self.Dim,1),[1,size(self.Dim,1)]),...
                                            [size(LabelCharges,1),1]),[size(self.Dim,1)*size(LabelCharges,1),1]);
                                        
                                        
                                        LabelCharges = [repmat(LabelCharges,[size(self.Dim,1),1]), reshape(repmat(reshape(1:size(self.Dim,1),[1,size(self.Dim,1)]),...
                                            [size(LabelCharges,1),1]),[size(self.Dim,1)*size(LabelCharges,1),1]), reshape(repmat(reshape(1:size(self.Dim,1),[1,size(self.Dim,1)]),...
                                            [size(LabelCharges,1),1]),[size(self.Dim,1)*size(LabelCharges,1),1])];
                                        
                                        
                                        if ExternalExtraChargeDirection ~= +1
                                            LabelCharges(:,end-1) = self.InverseIrrep(LabelCharges(:,end-1));
                                        end
                                        if InternalExtraChargeDirection ~= +1
                                            LabelCharges(:,end) = self.InverseIrrep(LabelCharges(:,end));
                                        end
                                        
                                        
                                        InputCharges = [InputCharges, ExternalExtraLeg,InternalExtraLeg];
                                        OutputCharges = [ExternalExtraLeg,InternalExtraLeg];
                                        
                                        Entries = repmat(Entries,[size(self.Dim,1),1]);
                                        
                                        %if self.IsMultiplicities
                                            LabelMults = [repmat(LabelMults,[size(self.Dim,1),1]), ones([size(self.Dim,1)*size(LabelMults,1),1])];
                                            
                                            %will only have one extra mult
                                            %value, as this must be the
                                            %trivial one
                                            
                                            MultiplicitiesExtra = WordDescriptions{kk}{ll,2}(4,1);
                                            InputMult = [InputMult, MultiplicitiesExtra];
                                        %end
                                        
                                        LabelChargesOutput = LabelCharges(:,[end-1,end]);
                                        OutputMax = self.MaxCharges.^(size(LabelChargesOutput,2));
                                        OutputSizes = self.MaxCharges.^(0:(size(LabelChargesOutput,2)-1));
                                        
                                    end
                                else
                                    LabelChargesOutput = LabelCharges(:,2*NumberChargesCombined);
                                    OutputMax = self.MaxCharges.^(size(LabelChargesOutput,2));
                                    OutputSizes = self.MaxCharges.^(0:(size(LabelChargesOutput,2)-1));
                                end
                                
                                OutputMult = zeros([1,0]);
                                %if ~self.IsMultiplicities
                                %    LabelMults = zeros([size(LabelCharges,1),0]);
                                %    InputMult = zeros([1,0]);
                                %end
                                
                                %now get rid of double ups:
                                
                                [~,Keep,~] = unique(InputCharges);
                                Keep = sort(Keep);
                                InputCharges = InputCharges(Keep);
                                LabelCharges = LabelCharges(:,Keep);
                                
                                %if self.IsMultiplicities
                                    [~,Keep,~] = unique(InputMult);
                                    Keep = sort(Keep);
                                    InputMult = InputMult(Keep);
                                    LabelMults = LabelMults(:,Keep);
                                %end
                                
                                
                                
                                LabelChargesInput = [LabelCharges, LabelMults];
                                
                                InputMax = self.MaxCharges.^(size(LabelCharges,2))*self.MaxMultiplicities.^(size(LabelMults,2));
                                InputSizes = [self.MaxCharges.^(0:(size(LabelCharges,2)-1)), self.MaxCharges.^(size(LabelCharges,2))*self.MaxMultiplicities.^(0:(size(LabelMults,2)-1))];
                                if isempty(InputCharges) && isempty(OutputCharges)
                                    LeftB = 1;
                                    RightB = 1;
                                    %OutputValues = 1;
                                    %InputValues = 1;
                                elseif isempty(OutputCharges)
                                    
                                    if ~isempty(InputCharges)
                                        LeftB = LabelChargesInput;
                                        %InputValues = sum((LabelChargesInput-1).*repmat(InputSizes,[size(LabelChargesInput,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(Labels,1),1]);
                                        %InputValues = ones([size(Labels,1),1]);
                                    end
                                    RightB = ones([size(LabelChargesInput,1),1]);
                                    %OutputValues = ones([size(LabelChargesInput,1),1]);
                                    
                                else %we know that OutputCharges is not empty
                                    
                                    RightB = LabelChargesOutput;
                                    %OutputValues = sum((LabelChargesOutput-1).*repmat(OutputSizes,[size(LabelChargesOutput,1),1]),2)+1;
                                    
                                    if ~isempty(InputCharges)
                                        LeftB = LabelChargesInput;
                                        %InputValues = sum((LabelChargesInput-1).*repmat(InputSizes,[size(LabelChargesInput,1),1]),2)+1;
                                    else
                                        LeftB = ones([size(LabelChargesOutput,1),1]);
                                        %InputValues = ones([size(LabelChargesOutput,1),1]);
                                    end
                                    
                                    
                                end
                                
                                EntriesB = Entries;
                                %Moves = sparse(InputValues,OutputValues,Entries, InputMax,OutputMax);
                                
                            otherwise
                                error('Error: An Invalid Word was used')
                                
                                
                                
                        end
                        
                        if ~FlagSkip %this allows us to ignore F-moves or Braiding if not used at all
                            InputRightCharges = InitialExtNumber + InitialIntNumber + InitialMultNumber + ...
                                [abs(InputCharges) + (InputCharges>0).*FinalExtNumber, FinalExtNumber + FinalIntNumber + InputMult]; 
                            %this is locations in lables that we pull to the right
                            InputLeftCharges = 1:(InitialExtNumber+InitialIntNumber+InitialMultNumber+FinalExtNumber+FinalIntNumber+FinalMultNumber);
                            InputLeftCharges(InputRightCharges) = [];
                            
                            
                            FinalExtNumber = FinalExtNumber - sum(InputCharges<0)+sum(OutputCharges<0);
                            FinalIntNumber = FinalIntNumber - sum(InputCharges>0)+sum(OutputCharges>0);
                            FinalMultNumber = FinalMultNumber - numel(InputMult)+numel(OutputMult);
                            
                            
                            %need to adjust to account for the lost legs
                            %which were subtracted:
                            OutputChargesOld = OutputCharges;
                            OutputMultOld = OutputMult;
                            for bb = 1:numel(OutputCharges)
                                if OutputCharges(bb)>0
                                    OutputCharges(bb) = OutputCharges(bb) + sum(OutputChargesOld(OutputChargesOld>0)<OutputCharges(bb))-sum(InputCharges(InputCharges>0)<OutputCharges(bb));
                                else %OutputCharges(bb)<0
                                    OutputCharges(bb) = OutputCharges(bb) - sum(OutputChargesOld(OutputChargesOld<0)>OutputCharges(bb))+sum(InputCharges(InputCharges<0)>OutputCharges(bb));
                                end
                            end
                            
                            
                            
                            OutputRightCharges = InitialExtNumber + InitialIntNumber + InitialMultNumber + [abs(OutputCharges) + (OutputCharges>0).*FinalExtNumber, FinalExtNumber + FinalIntNumber + OutputMult]; 
                            
                            
                            %this is locations in lables that we pull to the right
                            OutputLeftCharges = 1:(InitialExtNumber+InitialIntNumber+InitialMultNumber+FinalExtNumber+FinalIntNumber+FinalMultNumber);
                            OutputLeftCharges(OutputRightCharges) = [];
                            
                            %work out how to convert back to original
                            %charge list
                            [~,OutputPermute] = sort([OutputLeftCharges,OutputRightCharges], 'ascend');
                            
                            
                            
                            
                            
                            %Generate the sparse matrix S (note we need to
                            %use the max number to make sure the sizes are
                            %the same.
                            if ~isempty(InputLeftCharges)
                                LeftIndicies = Labels(:,InputLeftCharges);
                            else
                                LeftIndicies = ones([size(Labels,1),1]);
                            end
                            
                            if ~isempty(InputRightCharges)
                                RightIndicies = Labels(:,InputRightCharges);
                            else
                                RightIndicies = ones([size(Labels,1),1]);
                            end
                            
                            [LeftIndiciesOldLabels, ~,LeftIndiciesNum] = unique(LeftIndicies,'rows');
                            [RightIndiciesNewLabels, ~,RightBNum] = unique(RightB,'rows');
                            
                            
                            [~,~,MidNum] = unique([RightIndicies;LeftB],'rows');
                            RightIndiciesNum = MidNum(1:size(RightIndicies,1));
                            LeftBNum = MidNum((1+size(RightIndicies,1)):end);
                            
                            %multiply by the move sparse matrix  
                            SMatrix = sparse(LeftIndiciesNum, RightIndiciesNum, EntriesS, max(LeftIndiciesNum), max(MidNum));
                            BMatrix = sparse(LeftBNum,RightBNum,EntriesB, max(MidNum), max(RightBNum));
                            
                            SMatrix = SMatrix*BMatrix;
                            
                            [LeftIndiciesNum, RightIndiciesNum, EntriesS] = find(SMatrix);
                            
                            LeftIndicies = LeftIndiciesOldLabels(LeftIndiciesNum,:);
                            RightIndicies = RightIndiciesNewLabels(RightIndiciesNum,:);
                            
                            
                            %compute the new labels
                            if ~isempty(OutputLeftCharges)
                                if ~isempty(OutputRightCharges)
                                    Labels = [LeftIndicies,RightIndicies];
                                else
                                    Labels = LeftIndicies;
                                end
                            else
                                error('Affirmation Error: There should always be outputLeftCharges')
                            end
                            
                            %Labels = [LeftModIndicies,RightModIndicies];
                            Labels = Labels(:, OutputPermute);
                        end
                        
                    elseif isequal(WordDescriptions{kk}{ll,1}, 'Swap')
                        
                        
                        %the Permute will only be for relabeling 
                        
                        LabelsPermute = WordDescriptions{kk}{ll,2}(1,:);
                        LabelsPermuteMult = WordDescriptions{kk}{ll,2}(3,1:FinalMultNumber);
                        InvertLabel = WordDescriptions{kk}{ll,2}(2,:);
                        
                        LabelsPermute = abs(LabelsPermute)+(LabelsPermute>0).*FinalExtNumber;
                        MultPermute = LabelsPermuteMult+FinalIntNumber+FinalExtNumber;
                        
                        TrivLabels = isnan(LabelsPermute);
                        TrivMult =  isnan(MultPermute);
                        
                        LabelsWorking = Labels(:,(InitialExtNumber+InitialIntNumber+InitialMultNumber+1):end);
                        
                        LabelsWorking(:,[~TrivLabels,~TrivMult]) = LabelsWorking(:,[LabelsPermute(~TrivLabels),MultPermute(~TrivMult)]);
                        if any(TrivMult)
                            LabelsWorking(:,[TrivLabels,false(size(TrivMult))]) = self.TrivialIrrep;
                        end
                        if any(TrivLabels)
                            LabelsWorking(:,[false(size(TrivLabels)),TrivMult]) = 1;
                        end
                        
                        for nn = 1:size(InvertLabel,2)
                            if ~TrivLabels(nn)&&InvertLabel(nn)
                                LabelsWorking(:,nn) = self.InverseIrrep(LabelsWorking(:,nn));
                            end
                        end
                        
                        
                        Labels = [Labels(:,1:(InitialExtNumber + InitialIntNumber + InitialMultNumber)),LabelsWorking];
                        
                    elseif isequal(WordDescriptions{kk}{ll,1}, 'ForceSame')
                        Details = WordDescriptions{kk}{ll,2};
                        
                        DetailsLabels = Details(1:2,:);
                        
                        TrivialList = [DetailsLabels(1,~isnan(DetailsLabels(1,:))&isnan(DetailsLabels(2,:))), DetailsLabels(2,isnan(DetailsLabels(1,:))&~isnan(DetailsLabels(2,:)))];
                        
                        TrivialListExternal = TrivialList(TrivialList<0);
                        TrivialListInternal = TrivialList(TrivialList>0);
                        
                        Details = Details(:,any(isnan(DetailsLabels),1));
                        ChargeSame = Details(3,:).*Details(4,:);
                        Details = Details(1:2,:);
                        
                        ExternalDetailsKeep = all(Details<0,1);
                        
                        ExternalDetails = Details(:,ExternalDetailsKeep);
                        ExternalChargeSame = ChargeSame(:,ExternalDetailsKeep);
                        InternalDetails = Details(:,~ExternalDetailsKeep);
                        InternalChargeSame = ChargeSame(:,~ExternalDetailsKeep);
                        
                        TrivialCharges = [abs(TrivialListInternal), InitialExtNumber + InitialIntNumber + InitialMultNumber +[abs(TrivialListExternal), (TrivialListInternal)+FinalExtNumber]];
                        if ~isempty(TrivialCharges)
                            TrivialKeep = all(Labels(:,TrivialCharges)==self.TrivialIrrep,2);
                            Labels = Labels(TrivialKeep,:);
                            EntriesS = EntriesS(TrivialKeep);
                        end
                        
                        DetailsCharges1  = [abs(ExternalDetails(1,:)), InitialExtNumber + InitialIntNumber + InitialMultNumber +[abs(ExternalDetails(1,:)), (InternalDetails(1,:))+FinalExtNumber]];
                        DetailsCharges2  = [abs(ExternalDetails(2,:)), InitialExtNumber + InitialIntNumber + InitialMultNumber +[abs(ExternalDetails(2,:)), (InternalDetails(2,:))+FinalExtNumber]];
                        FullChargesSame = [ExternalChargeSame,ExternalChargeSame,InternalChargeSame];
                        
                        
                        if ~isempty(DetailsCharges1)
                            LabelsForSame1 = Labels(:,DetailsCharges1);
                            LabelsForSame2 = Labels(:,DetailsCharges2);
                            for uu = 1:numel(FullChargesSame)
                                if FullChargesSame(uu) == -1;
                                    LabelsForSame2(:,uu) = self.InverseIrrep(LabelsForSame2(:,uu));
                                end
                            end
                            
                            SameKeep = all(LabelsForSame1==LabelsForSame2,2);
                            Labels = Labels(SameKeep,:);
                            EntriesS = EntriesS(SameKeep);
                        end
                        
                    end
                    
                    if ll == 12
                        1;
                    end
                    %[ll,size(EntriesS)]
                end
                
                
                LabelsIn = Labels(:,1:(InitialExtNumber+InitialIntNumber+InitialMultNumber));
                LabelsOut = Labels(:,(InitialExtNumber+InitialIntNumber+InitialMultNumber+1):end);
                [~,LocationIn,ListIn] = unique([LabelsInAll;LabelsIn],'rows','first');
                %[~,Index] = sort(LocationIn,'ascend');
                ListIn = LocationIn(ListIn((size(LabelsInAll,1)+1):end));
                %[~,ListIn] = sort(ListIn,'ascend')
                [LabelsOut,~,ListOut] = unique(LabelsOut,'rows');
                
                OutputConvertNumbers{kk} = {LabelsInAll,LabelsOut,sparse(ListIn,ListOut,EntriesS,size(LabelsInAll,1),size(LabelsOut,1))};
                
            end
            
            if FlagOutputSingle
                OutputConvertNumbers = OutputConvertNumbers{1};
            end
            
        end
        
        
    end
    
    %SymCon save/load functions
    methods
        
        %checks if things have been saved or are valid:
        
        function BoolVal = CheckSaved(self, FlagNumber)
            if length(self.FlagSaved) < FlagNumber
                self.FlagSaved(FlagNumber) = false;
                self.FlagSavedMid(FlagNumber) = false;
                self.FlagSavedFine(FlagNumber) = false;
                self.SavedInfoDetails{FlagNumber,size(self.SavedInfoDetails,2)} = [];
                self.SavedCheckDetails{FlagNumber,size(self.SavedCheckDetails,2)} = [];
                self.SavedInfoMidDetails{FlagNumber,size(self.SavedInfoMidDetails,2)} = [];
                self.SavedCheckMidDetails{FlagNumber,size(self.SavedCheckMidDetails,2)} = [];
                self.SavedInfoFineDetails{FlagNumber,size(self.SavedInfoFineDetails,2)} = [];
                self.SavedCheckFineDetails{FlagNumber,size(self.SavedCheckFineDetails,2)} = [];
            end
            
            BoolVal = self.FlagSaved(FlagNumber);
        end
        
        function BoolVal = CheckSavedMid(self, FlagNumber)
            if length(self.FlagSaved) < FlagNumber
                self.FlagSaved(FlagNumber) = false;
                self.FlagSavedMid(FlagNumber) = false;
                self.FlagSavedFine(FlagNumber) = false;
                self.SavedInfoDetails{FlagNumber,size(self.SavedInfoDetails,2)} = [];
                self.SavedCheckDetails{FlagNumber,size(self.SavedCheckDetails,2)} = [];
                self.SavedInfoMidDetails{FlagNumber,size(self.SavedInfoMidDetails,2)} = [];
                self.SavedCheckMidDetails{FlagNumber,size(self.SavedCheckMidDetails,2)} = [];
                self.SavedInfoFineDetails{FlagNumber,size(self.SavedInfoFineDetails,2)} = [];
                self.SavedCheckFineDetails{FlagNumber,size(self.SavedCheckFineDetails,2)} = [];
            end
            
            BoolVal = self.FlagSavedMid(FlagNumber);
        end
        
        function BoolVal = CheckSavedFine(self, FlagNumber)
            if length(self.FlagSaved) < FlagNumber
                self.FlagSaved(FlagNumber) = false;
                self.FlagSavedMid(FlagNumber) = false;
                self.FlagSavedFine(FlagNumber) = false;
                self.SavedInfoDetails{FlagNumber,size(self.SavedInfoDetails,2)} = [];
                self.SavedCheckDetails{FlagNumber,size(self.SavedCheckDetails,2)} = [];
                self.SavedInfoMidDetails{FlagNumber,size(self.SavedInfoMidDetails,2)} = [];
                self.SavedCheckMidDetails{FlagNumber,size(self.SavedCheckMidDetails,2)} = [];
                self.SavedInfoFineDetails{FlagNumber,size(self.SavedInfoFineDetails,2)} = [];
                self.SavedCheckFineDetails{FlagNumber,size(self.SavedCheckFineDetails,2)} = [];
            end
            
            BoolVal = self.FlagSavedFine(FlagNumber);
        end
        
        function BoolVal = CheckValidLoad(self, FlagNumber,  Legs, ContractionOrder, LegOrder, ChargeSideTensorDetails,StructureTensorDetails,ChargeDimensionsTensorDetails)
            
            if length(self.FlagSaved) < FlagNumber
                self.FlagSaved(FlagNumber) = false;
                self.FlagSavedMid(FlagNumber) = false;
                self.FlagSavedFine(FlagNumber) = false;
                self.SavedInfoDetails{FlagNumber,size(self.SavedInfoDetails,2)} = [];
                self.SavedCheckDetails{FlagNumber,size(self.SavedCheckDetails,2)} = [];
                self.SavedInfoMidDetails{FlagNumber,size(self.SavedInfoMidDetails,2)} = [];
                self.SavedCheckMidDetails{FlagNumber,size(self.SavedCheckMidDetails,2)} = [];
                self.SavedInfoFineDetails{FlagNumber,size(self.SavedInfoFineDetails,2)} = [];
                self.SavedCheckFineDetails{FlagNumber,size(self.SavedCheckFineDetails,2)} = [];
            end
            
            BoolVal = self.FlagSaved(FlagNumber);
            
            BoolVal = BoolVal && (isequal(self.SavedCheckDetails{FlagNumber,1}, Legs) && isequal(self.SavedCheckDetails{FlagNumber,2}, ContractionOrder) &&...
                isequal(self.SavedCheckDetails{FlagNumber,3}, LegOrder)) && isequal(self.SavedCheckDetails{FlagNumber,4}, ChargeSideTensorDetails)&&...
                isequal(self.SavedCheckDetails{FlagNumber,5},StructureTensorDetails) && isequal(self.SavedCheckDetails{FlagNumber,6},ChargeDimensionsTensorDetails);
        end
        
        function BoolVal = CheckValidLoadMid(self, FlagNumber,  ExternalChargeTensorDetails,InternalChargeTensorDetails,MultiplicitiesInternalTensorDetails)
            BoolVal = self.FlagSavedMid(FlagNumber);
            BoolVal = BoolVal && (isequal(self.SavedCheckMidDetails{FlagNumber,1}, ExternalChargeTensorDetails) && isequal(self.SavedCheckMidDetails{FlagNumber,2}, InternalChargeTensorDetails) &&...
                isequal(self.SavedCheckMidDetails{FlagNumber,3}, MultiplicitiesInternalTensorDetails));
        end
        
        function BoolVal = CheckValidLoadFine(self, FlagNumber,  ChargeLegDimensionsTensorDetails)
            BoolVal = self.FlagSavedFine(FlagNumber);
            BoolVal = BoolVal && (isequal(self.SavedCheckFineDetails{FlagNumber,1}, ChargeLegDimensionsTensorDetails));
        end
        
        
        
        
        
        
        %clearing previously saved data:
        
        function ClearSaved(self,FlagNumber)
            if nargin<2
                self.FlagSaved = false([0,1]);
                self.SavedInfoDetails = cell([0,size(self.SavedInfoDetails,2)]);
                self.SavedCheckDetails = cell([0,size(self.SavedCheckDetails,2)]);
        
                self.FlagSavedMid = false([0,1]);
                self.SavedCheckMidDetails = cell([0,size(self.SavedCheckMidDetails,2)]);
                self.SavedInfoMidDetails = cell([0,size(self.SavedInfoMidDetails,2)]);
        
                self.FlagSavedFine = false([0,1]);
                self.SavedCheckFineDetails = cell([0,size(self.SavedCheckFineDetails,2)]);
                self.SavedInfoFineDetails = cell([0,size(self.SavedInfoFineDetails,2)]);
            else
                FlagNumber = unique(FlagNumber(:));
                FlagNumber = FlagNumber(FlagNumber<numel(self.FlagSaved));
                
                if ~isempty(FlagNumber)
                    self.FlagSaved(FlagNumber) = false([numel(FlagNumber),1]);
                    self.SavedInfoDetails(FlagNumber,:) = cell([numel(FlagNumber),size(self.SavedInfoDetails,2)]);
                    self.SavedCheckDetails(FlagNumber,:) = cell([numel(FlagNumber),size(self.SavedCheckDetails,2)]);
                end
                
                FlagNumber = FlagNumber(FlagNumber<numel(self.FlagSavedMid));
                if ~isempty(FlagNumber)
                    self.FlagSavedMid(FlagNumber) = false([numel(FlagNumber),1]);
                    self.SavedCheckMidDetails(FlagNumber,:) = cell([numel(FlagNumber),size(self.SavedCheckMidDetails,2)]);
                    self.SavedInfoMidDetails(FlagNumber,:) = cell([numel(FlagNumber),size(self.SavedInfoMidDetails,2)]);
                end
                
                
                FlagNumber = FlagNumber(FlagNumber<numel(self.FlagSavedFine));
                if ~isempty(FlagNumber)
                    self.FlagSavedFine(FlagNumber) = false([numel(FlagNumber),1]);
                    self.SavedCheckFineDetails(FlagNumber,:) = cell([numel(FlagNumber),size(self.SavedCheckFineDetails,2)]);
                    self.SavedInfoFineDetails(FlagNumber,:) = cell([numel(FlagNumber),size(self.SavedInfoFineDetails,2)]);
                end
            end
        end
        
        
        
        
        
        %save/load for SymCon:
        
        function SaveCurrentLoad(self, FlagNumber, CheckDetailsToSave, InfoDetailsToSave)
            
            %{Legs, ContractionOrder, LegOrder} = CheckDetailsToSave(:)';
            %{ActionTensorsLocations, ActualLegsTotal,ActionTypesTotal, LegsContractedTotal, KeepUnique, OldTrueLocations, PermuteFinal} = InfoDetailsToSave(:)';
            
            %HERAFix:  need to put checks in:
            self.FlagSaved(FlagNumber) = true;
            self.SavedCheckDetails(FlagNumber,:) = CheckDetailsToSave(:)';
            self.SavedInfoDetails(FlagNumber,:) = InfoDetailsToSave(:)';
            
        end
        
        function SaveCurrentLoadMid(self, FlagNumber, CheckMidDetailsToSave, InfoMidDetailsToSave)
            
            %{Legs, ContractionOrder, LegOrder} = CheckMidDetailsToSave(:)';
            %{ActionTensorsLocations, ActualLegsTotal,ActionTypesTotal, LegsContractedTotal, KeepUnique, OldTrueLocations, PermuteFinal} = InfoMidDetailsToSave(:)';
            
            %HERAFix:  need to put checks in:
            self.FlagSavedMid(FlagNumber) = true;
            self.SavedCheckMidDetails(FlagNumber,:) = CheckMidDetailsToSave(:)';
            self.SavedInfoMidDetails(FlagNumber,:) = InfoMidDetailsToSave(:)';
            
        end
        
        function SaveCurrentLoadFine(self, FlagNumber, CheckFineDetailsToSave, InfoFineDetailsToSave)
            
            %{Legs, ContractionOrder, LegOrder} = CheckFineDetailsToSave(:)';
            
            %HERAFix:  need to put checks in:
            self.FlagSaved(FlagNumber) = true;
            self.SavedCheckFineDetails(FlagNumber,:) = CheckFineDetailsToSave(:)';
            %self.SavedInfoDetails(FlagNumber,:) = InfoFineDetailsToSave(:)';
            
        end
        
        function SaveCurrentLoadPlanar(self, FlagNumber, CheckPlanarDetailsToSave, InfoPlanarDetailsToSave)
            
            %{Legs, ContractionOrder, LegOrder} = CheckFineDetailsToSave(:)';
            
            %HERAFix:  need to put checks in:
            self.FlagSavedPlanar(FlagNumber) = true;
            %self.SavedCheckPlanarDetails(FlagNumber,:) = CheckPlanarDetailsToSave(:)';
            self.SavedInfoPlanarDetails(FlagNumber,:) = InfoPlanarDetailsToSave(:)';
            
        end
        
        
        
        
        
        
        
        function [ActionTensorsLocations, ActionTypeTotal, ActionTensorsPermute, ActionNumberContractedLegs,...
                 OldTrueLocations, OldTruePermute, OldTrueChargeSides, OldTrueChargeSidesInt,...
                 OldTrueStructure,OldTrueChargeDirections,ActionTotalTensors,DropTensor,OldTrue,OldTrueBraiding,OldTrueBraidingDirection] = LoadCurrentData(self, FlagNumber)
            ActionTensorsLocations = self.SavedInfoDetails{FlagNumber,1};
            ActionTypeTotal = self.SavedInfoDetails{FlagNumber,2};
            ActionTensorsPermute = self.SavedInfoDetails{FlagNumber,3};
            ActionNumberContractedLegs = self.SavedInfoDetails{FlagNumber,4};
            OldTrueLocations = self.SavedInfoDetails{FlagNumber,5};
            OldTruePermute = self.SavedInfoDetails{FlagNumber,6};
            OldTrueChargeSides = self.SavedInfoDetails{FlagNumber,7};
            OldTrueChargeSidesInt = self.SavedInfoDetails{FlagNumber,8};
            OldTrueStructure = self.SavedInfoDetails{FlagNumber,9};
            OldTrueChargeDirections = self.SavedInfoDetails{FlagNumber,10};
            
            ActionTotalTensors = self.SavedInfoDetails{FlagNumber,11};
            DropTensor = self.SavedInfoDetails{FlagNumber,12};
            OldTrue = self.SavedInfoDetails{FlagNumber,13};
            OldTrueBraiding = self.SavedInfoDetails{FlagNumber,14};
            OldTrueBraidingDirection = self.SavedInfoDetails{FlagNumber,15};
        end
        
        function [UniqueLegsNeg, FlagUserDidNotChooseLegOrder] = LoadDataForPlanarWithoutBroad(self, FlagNumber)
            UniqueLegsNeg = self.SavedInfoDetails{FlagNumber,16};
            FlagUserDidNotChooseLegOrder = self.SavedInfoDetails{FlagNumber,17};
        end
        
        function [TotalTensors, ActualLegsTotal, LegOrder,ActionTensorsStructure,ActionTensorsChargeDirectionsExternalOriginal] = LoadDataForMidWithoutBroad(self, FlagNumber)
            TotalTensors = self.SavedInfoDetails{FlagNumber,18};
            ActualLegsTotal = self.SavedInfoDetails{FlagNumber,19};
            LegOrder = self.SavedInfoDetails{FlagNumber,20};
            ActionTensorsStructure = self.SavedInfoDetails{FlagNumber,21};
            ActionTensorsChargeDirectionsExternalOriginal = self.SavedInfoDetails{FlagNumber,22};
        end

        function [UniqueLegsPos] = LoadDataForFineWithoutBroad(self, FlagNumber)
            UniqueLegsPos = self.SavedInfoDetails{FlagNumber,23};
        end
        
        
        
        
        
        
        function [C_MultiplicitiesInternal, C_ChargeLabelsInternal, C_ChargeLabelsExternal, ActionChargesList, C_List,ActionMaxNumTensorParts,...
                OldTrueChargesExternal,OldTrueChargesInternal,OldTrueMultiplicitiesInternal] = LoadCurrentMidData(self, FlagNumber)
            C_MultiplicitiesInternal = self.SavedInfoMidDetails{FlagNumber,1};
            C_ChargeLabelsInternal = self.SavedInfoMidDetails{FlagNumber,2};
            C_ChargeLabelsExternal = self.SavedInfoMidDetails{FlagNumber,3};
            ActionChargesList = self.SavedInfoMidDetails{FlagNumber,4};
            C_List = self.SavedInfoMidDetails{FlagNumber,5};
            ActionMaxNumTensorParts = self.SavedInfoMidDetails{FlagNumber,6};
            OldTrueChargesExternal = self.SavedInfoMidDetails{FlagNumber,7};
            OldTrueChargesInternal = self.SavedInfoMidDetails{FlagNumber,8};
            OldTrueMultiplicitiesInternal = self.SavedInfoMidDetails{FlagNumber,9};
        end
        
        
        function [OldTrueConvertWords, ActionTensorsConvertWords, ActionTensorsCombineCharges,ActionTensorsChargeDirectionsExternalOriginal] = LoadDataForMidWithoutPlanar(self, FlagNumber)
            OldTrueConvertWords = self.SavedInfoPlanarDetails{FlagNumber,1};
            ActionTensorsConvertWords = self.SavedInfoPlanarDetails{FlagNumber,2};
            ActionTensorsCombineCharges = self.SavedInfoPlanarDetails{FlagNumber,3};
            ActionTensorsChargeDirectionsExternalOriginal = self.SavedInfoPlanarDetails{FlagNumber,8};
        end

        
        function [OldTrueConvertNumbers, ActionTensorsConvertWords, ActionTensorsCombineCharges, OldTrueListMultNumbers, ActionTensorsListMultNumbers,...
                OldTrueList,ActionChargesList,ActionTensorsChargeDirectionsExternalOriginal] = LoadCurrentPlanarMidData(self, FlagNumber)
            OldTrueConvertNumbers = self.SavedInfoPlanarDetails{FlagNumber,1};
            ActionTensorsConvertWords = self.SavedInfoPlanarDetails{FlagNumber,2};
            ActionTensorsCombineCharges = self.SavedInfoPlanarDetails{FlagNumber,3};
            OldTrueListMultNumbers = self.SavedInfoPlanarDetails{FlagNumber,4};
            ActionTensorsListMultNumbers = self.SavedInfoPlanarDetails{FlagNumber,5};
            OldTrueList = self.SavedInfoPlanarDetails{FlagNumber,6};
            ActionChargesList = self.SavedInfoPlanarDetails{FlagNumber,7};
            ActionTensorsChargeDirectionsExternalOriginal = self.SavedInfoPlanarDetails{FlagNumber,8};
        end

        

        
        
        %External save/load functions for between input parameters:
        
        function SaveExternal(self, ExternalFileName)
            
            RealName = self.Name;
            
            RealFlagSaved = self.FlagSaved;
            RealCheckDetails = self.SavedCheckDetails;
            RealInfoDetails = self.SavedInfoDetails;
            
            MidFlagSaved = self.FlagSavedMid;
            MidCheckDetails = self.SavedCheckMidDetails;
            MidInfoDetails = self.SavedInfoMidDetails;
            
            FineFlagSaved = self.FlagSavedFine;
            FineCheckDetails = self.SavedCheckFineDetails;
            FineInfoDetails = self.SavedInfoFineDetails;
            
            PlanarFlagSaved = self.FlagSavedPlanar;
            PlanarCheckDetails = self.SavedCheckPlanarDetails;
            PlanarInfoDetails = self.SavedInfoPlanarDetails;
            
            MidPlanarFlagSaved = self.FlagSavedMidPlanar;
            
            save(ExternalFileName, 'RealName','RealFlagSaved', 'RealCheckDetails','RealInfoDetails',...
                'MidFlagSaved', 'MidCheckDetails', 'MidInfoDetails','FineFlagSaved','FineCheckDetails','FineInfoDetails',...
                'PlanarFlagSaved','PlanarCheckDetails','PlanarInfoDetails','MidPlanarFlagSaved', '-v7.3');
        end
        
        function LoadExternal(self, ExternalFileName)
            load(ExternalFileName, 'RealName','RealFlagSaved', 'RealCheckDetails','RealInfoDetails',...
                'MidFlagSaved', 'MidCheckDetails', 'MidInfoDetails','FineFlagSaved','FineCheckDetails','FineInfoDetails',...
                'PlanarFlagSaved','PlanarCheckDetails','PlanarInfoDetails','MidPlanarFlagSaved');
            
            if isequal(self.Name, RealName)
                
                %HERAFix: need to add checks in
                
                self.FlagSaved = RealFlagSaved;
                self.SavedCheckDetails = RealCheckDetails;
                self.SavedInfoDetails = RealInfoDetails;
                
                self.FlagSavedMid = MidFlagSaved;
                self.SavedCheckMidDetails = MidCheckDetails;
                self.SavedInfoMidDetails = MidInfoDetails;
                
                self.FlagSavedFine = FineFlagSaved;
                self.SavedCheckFineDetails = FineCheckDetails;
                self.SavedInfoFineDetails = FineInfoDetails;
                
                self.FlagSavedPlanar = PlanarFlagSaved;
                self.SavedCheckPlanarDetails = PlanarCheckDetails;
                self.SavedInfoPlanarDetails = PlanarInfoDetails;
                
                MidPlanarFlagSaved = self.FlagSavedMidPlanar;
            else
                error('Loading was done from a different Symmetry type')
            end
            
        end

    end
    
    %Symmetry generation
    methods
        
        function self = Symmetry(label)
            
            
            FlagGeneratedSymmetry = false;
            
            %work out if there are any times
            TimesLocations = strfind(label, ' X ');
            PlusLocations = strfind(label, ' + ');
            OpenBracketLocations = strfind(label, '[');
            CloseBracketLocations = strfind(label, ']');
            
            %currently just using times:
            
            if ~isempty(TimesLocations)
                
                Alabel = label(1:(TimesLocations(1)-1));
                Blabel = label((TimesLocations(1)+3):end);
                
                ASym = Symmetry(Alabel);
                BSym = Symmetry(Blabel);
                
                %work out dimensions
                SizeDimA = size(ASym.Dim)-[0,2];
                SizeDimB = size(BSym.Dim)-[0,2];
                
                HumanNames = [reshape(repmat(reshape(ASym.Dim(:,3:end),[SizeDimA(1),1,SizeDimA(2)]),[SizeDimB(1),1]),[SizeDimA(1)*SizeDimB(1),SizeDimA(2)]),...
                    reshape(repmat(reshape(BSym.Dim(:,3:end),[1,SizeDimB(1),SizeDimB(2)]),[SizeDimA(1),1]),[SizeDimA(1)*SizeDimB(1),SizeDimB(2)])];
                SymDims = repmat(ASym.Dim(:,2),[SizeDimB(1),1]).*reshape(repmat(reshape(BSym.Dim(:,2),[1,SizeDimA(1),1]),[SizeDimB(1),1]),[SizeDimA(1)*SizeDimB(1),1]);
                
                MachineNamesA = repmat(ASym.Dim(:,1),[SizeDimB(1),1]);
                MachineNamesB = reshape(repmat(reshape(BSym.Dim(:,1),[1,SizeDimB(1),1]),[SizeDimA(1),1]),[SizeDimA(1)*SizeDimB(1),1]);
                
                self.Dim = [(1:(SizeDimA(1)*SizeDimB(1)))', SymDims, HumanNames];
                
                
                %work out Fusion2
                
                Fusion2OrigA = ASym.Fusion2;
                
                Fusion2ModA = repmat(ASym.Fusion2,[SizeDimB(1),SizeDimB(1),SizeDimB(1)]);
                Fusion2ModB = repmat(BSym.Fusion2,[SizeDimA(1),SizeDimA(1),SizeDimA(1)]);
                
                RelabelB = reshape(repmat(1:SizeDimB(1),[SizeDimA(1),1])+repmat((0:(SizeDimA(1)-1))',[1,SizeDimB(1)])*2,[SizeDimA(1)*SizeDimB(1),1]);
                
                Fusion2ModB = Fusion2ModB(RelabelB,RelabelB,RelabelB);
                
                self.Fusion2 = Fusion2ModA.*Fusion2ModB;
                
                %work out F-moves
                
                NumberStoredA = size(ASym.F{1});
                NumberStoredB = size(BSym.F{1});
                
                FOriginalA = reshape(repmat(reshape(ASym.F{1},[NumberStoredA(1),1,NumberStoredA(2)]),[NumberStoredB(1),1]),[NumberStoredA(1)*NumberStoredB(1),NumberStoredA(2)]);
                FOriginalB = reshape(repmat(reshape(BSym.F{1},[1,NumberStoredB(1),NumberStoredB(2)]),[NumberStoredA(1),1]),[NumberStoredA(1)*NumberStoredB(1),NumberStoredB(2)]);
                
                FMixed = FOriginalA(:,1:6)+(FOriginalB(:,1:6)-1)*SizeDimA(1);
                
                if NumberStoredA(2) == 10 && NumberStoredB(2) == 10
                    ReleventMultALocations = [FOriginalA(:,1)];%sum together FOriginal to pull out entry in Fusion2
                    
                    ReleventMultA = Fusion2OrigA(ReleventMultALocations);
                    
                    FMixed = [FMixed,FOriginalA(:,7:10)+(FOriginalB(:,7:10)-1).*ReleventMultA];
                    
                elseif NumberStoredA(2) == 10 
                    FMixed = [FMixed, FOriginalA(:,7:10)];
                    
                elseif NumberStoredB(2) == 10
                    FMixed = [FMixed, FOriginalB(:,7:10)];
                    
                end
                
                self.F{1} = FMixed;
                
                self.F{2} = repmat(ASym.F{2},[NumberStoredB(1),1]).*reshape(repmat(reshape(BSym.F{2},[1,NumberStoredB(1),1]),[NumberStoredA(1),1]),[NumberStoredA(1)*NumberStoredB(1),1]);
                self.F{3} = repmat(ASym.F{3},[NumberStoredB(1),1]).*reshape(repmat(reshape(BSym.F{3},[1,NumberStoredB(1),1]),[NumberStoredA(1),1]),[NumberStoredA(1)*NumberStoredB(1),1]);
                
                
                %work out R-moves
                
                NumberStoredA = size(ASym.R{1});
                NumberStoredB = size(BSym.R{1});
                
                ROriginalA = reshape(repmat(reshape(ASym.R{1},[NumberStoredA(1),1,NumberStoredA(2)]),[NumberStoredB(1),1]),[NumberStoredA(1)*NumberStoredB(1),NumberStoredA(2)]);
                ROriginalB = reshape(repmat(reshape(BSym.R{1},[1,NumberStoredB(1),NumberStoredB(2)]),[NumberStoredA(1),1]),[NumberStoredA(1)*NumberStoredB(1),NumberStoredB(2)]);
                
                RMixed = ROriginalA(:,1:3)+(ROriginalB(:,1:3)-1)*SizeDimA(1);
                
                if NumberStoredA(2) == 5 && NumberStoredB(2) == 5
                    ReleventMultALocations = [FOriginalA(:,1)];%sum together FOriginal to pull out entry in Fusion2
                    
                    ReleventMultA = Fusion2OrigA(ReleventMultALocations);
                    
                    RMixed = [RMixed,ROriginalA(:,4:5)+(ROriginalB(:,4:5)-1).*ReleventMultA];
                    
                elseif NumberStoredA(2) == 5 
                    RMixed = [RMixed, ROriginalA(:,4:5)];
                    
                elseif NumberStoredB(2) == 5
                    RMixed = [RMixed, ROriginalB(:,4:5)];
                    
                end
                
                self.R{1} = RMixed;
                
                self.R{2} = repmat(ASym.R{2},[NumberStoredB(1),1]).*reshape(repmat(reshape(BSym.R{2},[1,NumberStoredB(1),1]),[NumberStoredA(1),1]),[NumberStoredA(1)*NumberStoredB(1),1]);
                self.R{3} = repmat(ASym.R{3},[NumberStoredB(1),1]).*reshape(repmat(reshape(BSym.R{3},[1,NumberStoredB(1),1]),[NumberStoredA(1),1]),[NumberStoredA(1)*NumberStoredB(1),1]);
                
                self.FlagBraiding = ASym.FlagBraiding||BSym.FlagBraiding;
                
                self.Name = [ASym.Name, ' X ', BSym.Name];
                
                FlagGeneratedSymmetry = true;
                
            end
            
            
            
            if ~FlagGeneratedSymmetry
            %temporarly I'm only going to do none, Z2, Z3, U(1), SU(2) or Fibonacci
%            self.FlagSym = true;
            if length(label)>=5 && isequal(label(1:5),'Refl(')
                FlagGeneratedSymmetry = true;
                
                OpenBrackets = find(label(6:end)=='(');
                OpenBracketCount = zeros([1,numel(label)-5]);
                if ~isempty(OpenBrackets)
                    OpenBracketCount(OpenBrackets) = 1;
                end
                OpenBracketCount = cumsum(OpenBracketCount);
                
                CloseBrackets = find(label(6:end)==')');
                CloseBracketCount = zeros([1,numel(label)-5]);
                if ~isempty(CloseBrackets)
                    CloseBracketCount(CloseBrackets) = 1;
                end
                CloseBracketCount = cumsum(CloseBracketCount);
                
                
                
                SumsBrackets = ones([1,numel(label)-5])+OpenBracketCount-CloseBracketCount;
                Loc = find(SumsBrackets==0, 1,'first');
                if isempty(Loc)
                    Loc = numel(label)+1;
                end
                
                TempSym = Symmetry(label(6:(Loc-1+5)));
                
                WorkingFusion2 = TempSym.Fusion2;
                WorkingDim = TempSym.Dim;
                WorkingF = TempSym.F;
                WorkingR = TempSym.R;
                
                
                A = cat(1,cat(2,WorkingFusion2,zeros(size(WorkingFusion2))),cat(2,zeros(size(WorkingFusion2)),WorkingFusion2));
                B = cat(1,cat(2,zeros(size(WorkingFusion2)),WorkingFusion2),cat(2,WorkingFusion2,zeros(size(WorkingFusion2))));
                self.Fusion2 = cat(3,A,B);
                
                TempNames = cell([4,1]);
                List = [0,0,0;...
                        1,1,0;...
                        1,0,1;...
                        0,1,1]';
                    
                Counter = 0;
                for ll = List
                    Counter = Counter+1;
                    TempNames{Counter} = WorkingF{1};
                    
                    TempNames{Counter}(:,1:3) = TempNames{Counter}(:,1:3)+size(WorkingDim,1)*repmat(ll',[size(WorkingF{1},1),1]);
                end
                
                WorkingR{1} = cell2mat(TempNames);
                WorkingR{2} = repmat(WorkingR{2},[4,1]);
                WorkingR{3} = repmat(WorkingR{3},[4,1]);
                
                
                TempNames = cell([8,1]);
                List = [0,0,0,0;...
                        0,0,1,1;...
                        0,1,0,1;...
                        0,0,1,1;...
                        1,0,1,0;...
                        1,1,0,0;...
                        0,1,1,0;...
                        1,1,1,1]';
                
                List = [List;mod([List(2,:);List(2,:)]+[List(3,:);List(1,:)],2)];
                
                Counter = 0;
                for ll = List
                    Counter = Counter+1;
                    TempNames{Counter} = WorkingF{1};
                    
                    TempNames{Counter}(:,1:6) = TempNames{Counter}(:,1:6)+size(WorkingDim,1)*repmat(ll',[size(WorkingF{1},1),1]);
                end
                
                WorkingF{1} = cell2mat(TempNames);
                WorkingF{2} = repmat(WorkingF{2},[8,1]);
                WorkingF{3} = repmat(WorkingF{3},[8,1]);
                
                
                self.Dim = [[WorkingDim,ones([size(WorkingDim,1),1])]; [WorkingDim,-ones([size(WorkingDim,1),1])]];
                self.Dim(:,1) = 1:size(self.Dim,1);
                
                self.Name = label;
                self.R = WorkingR;
                self.F = WorkingF;
            end

            if length(label)>=4 && isequal(label(1:4),'none') %this is for basic stuff
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0];
                self.Name = 'none';
                self.TrivialIrrep = 1;
                self.InverseIrrep = 1;
                self.Fusion2 = 1;
                
                self.F = 1;
                self.R = 1;
                self.FlagSphericalSymmetry = true;
                self.FlagSym = false;
            end

            if length(label)>=4 && isequal(label(1:4),'aZ3_') %this is for non-finite groups
                FlagGeneratedSymmetry = true;
                
                N = str2num(label(5:end));
                
                if N~=round(N)
                    error('error: need an integer value for the anyonic Zn models')
                end
                
                if N < 0 || N>=3
                    error('error: We need the counter to be Zn valued for the anyonic Zn models')
                end
                self.Dim = [1,1,0;2,1,1;3,1,2];
                self.Name = ['aZ3_',num2str(N)];
                self.TrivialIrrep = [1];
                self.InverseIrrep = [1;2;3];
                self.Fusion2 = cat(3,[1,0,0;0,1,0;0,0,1],[0,0,1;1,0,0;0,1,0],[0,1,0;0,0,1;1,0,0]);
                
                
                self.F = zeros(size(self.Dim,1)*[1,1,1,1,1,1]);
                for aa = 1:size(self.Dim,1)
                    for bb = 1:size(self.Dim,1)
                        for cc = 1:size(self.Dim,1)
                            
                            ee = find(self.Fusion2(:,aa,bb),1);
                            ff = find(self.Fusion2(:,bb,cc),1);
                            dd = find(self.Fusion2(:,ee,cc),1);
                            
                            self.F(aa,bb,cc,dd,ee,ff) = 1;
                            
                        end
                    end
                end
                
                self.R = zeros(size(self.Dim,1)*[1,1,1]);
                for aa = 1:size(self.Dim,1)
                    for bb = 1:size(self.Dim,1)
                        
                        cc = find(self.Fusion2(:,aa,bb),1);
                        
                        self.R(aa,bb,cc) = exp(1i*2*pi*N/3*(aa-1)*(bb-1));
                        
                    end
                end
                
            end
            
            if isequal(label,'Z2') %this is standard
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;2,1,1];
                self.Name = 'Z2';
                self.TrivialIrrep = [1];
                self.InverseIrrep = [1;2];
                self.Fusion2 = cat(3,[1,0;0,1],[0,1;1,0]);
            end
            
            if isequal(label,'Z3') %this is for non-self dual
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;2,1,1;3,1,2];
                self.Name = 'Z3';
                self.TrivialIrrep = [1];
                self.InverseIrrep = [1;3;2];
                self.Fusion2 = cat(3,[1,0,0;0,1,0;0,0,1],[0,0,1;1,0,0;0,1,0],[0,1,0;0,0,1;1,0,0]);
            end
            
            if isequal(label,'Fib') %this is for anyonic (non-abelian)
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;2,(1+sqrt(5))/2,1];
                self.Name = 'Fib';
                self.Fusion2 = cat(3,[1,0;0,1],[0,1;1,1]);
                
                phi = (1+sqrt(5))/2; % the golden ratio
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:5
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessFLabels = zeros([size(self.Dim,1)^4,6]);
                AccessFValues = zeros([size(self.Dim,1)^4,1]);
                AccessFInvValues = zeros([size(self.Dim,1)^4,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    ee = Work(5);
                    if self.Fusion2(ee,aa,bb)>0 && self.Fusion2(dd,ee,cc)>0
                        for ff = 1:size(self.Dim,1)
                            if self.Fusion2(ff,bb,cc)>0 && self.Fusion2(dd,aa,ff)>0
                                if any([aa,bb,cc,dd] ~= 2)
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = 1;
                                    AccessFInvValues(Counter,:) = 1;
                                    Counter = Counter+1;
                                else
                                    if ee == 1 && ff == 1
                                        AccessFLabels(Counter+(0:3),:) = [repmat([aa,bb,cc,dd],[4,1]),[1;1;2;2],[1;2;1;2]];
                                        Temp = phi.^(-[2,1;1,2]/2).*[1,1,;1,-1];
                                        TempInv = Temp^(-1);
                                        AccessFValues(Counter+(0:3),:) = Temp(:);
                                        AccessFInvValues(Counter+(0:3),:) = TempInv(:);
                                        Counter = Counter+4;
                                    end
                                end
                            end
                        end
                    end
                end
                
                self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:3
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessRLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRValues = zeros([size(self.Dim,1)^3,1]);
                AccessRInvValues = zeros([size(self.Dim,1)^3,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    if self.Fusion2(cc,aa,bb)>0
                        if any([aa,bb]~=2)
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = 1;
                            AccessRInvValues(Counter,1) = 1;
                            Counter = Counter+1;
                        else
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = exp(-1i*4*pi/5)*(cc==1) + exp(1i*3*pi/5)*(cc==2);
                            AccessRInvValues(Counter,1) = exp(1i*4*pi/5)*(cc==1) + exp(-1i*3*pi/5)*(cc==2);
                            Counter = Counter+1;
                        end
                    end
                end
                self.R = {AccessRLabels(1:(Counter-1),:), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                self.FlagBraiding = true;
                
            end
            
            if isequal(label,'Fermion') %this is for anyonic (non-abelian)
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;2,(1+sqrt(5))/2,1];
                self.Name = 'Fermion';
                self.Fusion2 = cat(3,[1,0;0,1],[0,1;1,0]);
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:5
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessFLabels = zeros([size(self.Dim,1)^4,6]);
                AccessFValues = zeros([size(self.Dim,1)^4,1]);
                AccessFInvValues = zeros([size(self.Dim,1)^4,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    ee = Work(5);
                    if self.Fusion2(ee,aa,bb)>0 && self.Fusion2(dd,ee,cc)>0
                       for ff = 1:size(self.Dim,1)
                           if self.Fusion2(ff,bb,cc)>0 && self.Fusion2(dd,aa,ff)>0
                                if any([aa,bb,cc,dd] ~= 2)
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = 1;
                                    AccessFInvValues(Counter,:) = 1;
                                    Counter = Counter+1;
                                else
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = 1;
                                    AccessFInvValues(Counter,:) = 1;
                                    Counter = Counter+1;
                                end               
                            end
                        end
                    end
                end
                
                self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:3
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessRLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRValues = zeros([size(self.Dim,1)^3,1]);
                AccessRInvValues = zeros([size(self.Dim,1)^3,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    if self.Fusion2(cc,aa,bb)>0
                        if aa == 2 && bb == 2 && cc == 1
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = -1;
                            AccessRInvValues(Counter,1) = -1;
                            Counter = Counter+1;
                        else
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = 1;
                            AccessRInvValues(Counter,1) = 1;
                            Counter = Counter+1;
                        end
                    end
                end
                self.R = {AccessRLabels(1:(Counter-1),:), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                self.FlagBraiding = true;
                
            end
            
            
            if isequal(label,'Ising') %this is for anyonic (non-abelian)
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;2,sqrt(2),1;3,1,2];
                self.Name = 'Ising';
                self.Fusion2 = cat(3,[1,0,0;0,1,0;0,0,1],[0,1,0;1,0,1;0,1,0],[0,0,1;0,1,0;1,0,0]);
                
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:5
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessFLabels = zeros([size(self.Dim,1)^4,6]);
                AccessFValues = zeros([size(self.Dim,1)^4,1]);
                AccessFInvValues = zeros([size(self.Dim,1)^4,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    ee = Work(5);
                    if self.Fusion2(ee,aa,bb)>0 && self.Fusion2(dd,ee,cc)>0
                        for ff = 1:size(self.Dim,1)
                            if self.Fusion2(ff,bb,cc)>0 && self.Fusion2(dd,aa,ff)>0
                                if all([aa,bb,cc,dd] == [2,3,2,3])
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = -1;
                                    AccessFInvValues(Counter,:) = -1;
                                    Counter = Counter+1;
                                elseif all([aa,bb,cc,dd] == [3,2,3,2])
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = -1;
                                    AccessFInvValues(Counter,:) = -1;
                                    Counter = Counter+1;
                                elseif all([aa,bb,cc,dd] == [2,2,2,2])
                                    if ee == 1 && ff == 1
                                        AccessFLabels(Counter+(0:3),:) = [repmat([aa,bb,cc,dd],[4,1]),[1;1;3;3],[1;3;1;3]];
                                        Temp = sqrt(2)^(-1)*[1,1,;1,-1];
                                        TempInv = Temp^(-1);
                                        AccessFValues(Counter+(0:3),:) = Temp(:);
                                        AccessFInvValues(Counter+(0:3),:) = TempInv(:);
                                        Counter = Counter+4;
                                    end
                                else
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = 1;
                                    AccessFInvValues(Counter,:) = 1;
                                    Counter = Counter+1;
                                end
                            end
                        end
                    end
                end
                
                self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:3
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessRLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRValues = zeros([size(self.Dim,1)^3,1]);
                AccessRInvValues = zeros([size(self.Dim,1)^3,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    if self.Fusion2(cc,aa,bb)>0
                        if all([aa,bb,cc] == [2,2,1])
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = exp(-1i*pi/8);
                            AccessRInvValues(Counter,1) = exp(1i*pi/8);
                            Counter = Counter+1;
                        elseif all([aa,bb,cc] == [2,2,3])
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = exp(1i*3*pi/8);
                            AccessRInvValues(Counter,1) = exp(-1i*3*pi/8);
                            Counter = Counter+1;
                        elseif all([aa,bb,cc] == [3,3,1])
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = -1;
                            AccessRInvValues(Counter,1) = -1;
                            Counter = Counter+1;
                        elseif all([aa,bb,cc] == [2,3,2])
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = exp(-1i*pi/2);
                            AccessRInvValues(Counter,1) = exp(1i*pi/2);
                            Counter = Counter+1;
                        elseif all([aa,bb,cc] == [3,2,2])
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = exp(-1i*pi/2);
                            AccessRInvValues(Counter,1) = exp(1i*pi/2);
                            Counter = Counter+1;
                        else
                            AccessRLabels(Counter,1:3) = [aa,bb,cc];
                            AccessRValues(Counter,1) = 1;
                            AccessRInvValues(Counter,1) = 1;
                            Counter = Counter+1;
                        end
                    end
                end
                self.R = {AccessRLabels(1:(Counter-1),:), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                self.FlagBraiding = true;
                
            end
            
            
            if numel(label)>=5 && isequal(label(1:5),'SU2k_') %this is for anyonic (non-abelian)
                FlagGeneratedSymmetry = true;
                
                k = str2num(label(6:end));
                
                if ~IsInteger(k)
                    error('error: We need a positive number for k in SU(2)k')
                end
                
                if k <= 0
                    error('error: We need a positive number for k in SU(2)k')
                end
                
                q = exp(1i*2*pi()/(k+2));
                
                Dims = sin((1:(k+1)')*pi()/(k+2))./sin(pi()/(k+2));
                self.Dim = [(1:(k+1))',  Dims(:), (0:k)'./2];
                self.Name = ['SU(2)_k=',num2str(k)];
                self.TrivialIrrep = [1];
                self.InverseIrrep = (1:(k+1))';
                
                FusionPart = cell([1,k+1]);
                for jj2 = 1:(k+1)
                    Part = zeros([k+1,k+1]);
                    for jj1 = 1:(k+1)
                        
                        UseMin = abs(jj1-jj2)+1;
                        UseMax = min(jj1+jj2-1,2*k+2-(jj1+jj2-1));
                        
                        if UseMax>=UseMin
                            Part(UseMin:2:UseMax,jj1) = 1;
                        end
                    end
                    FusionPart{jj2} = Part;
                end
                
                self.Fusion2 = cat(3,FusionPart{:});
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:5
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessFLabels = zeros([size(self.Dim,1)^4,6]);
                AccessFValues = zeros([size(self.Dim,1)^4,1]);
                AccessFInvValues = zeros([size(self.Dim,1)^4,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    ee = Work(5);
                    if self.Fusion2(ee,aa,bb)>0 && self.Fusion2(dd,ee,cc)>0
                        for ff = 1:size(self.Dim,1)
                            if self.Fusion2(ff,bb,cc)>0 && self.Fusion2(dd,aa,ff)>0
                                jaa = (aa-1)/2;
                                jbb = (bb-1)/2;
                                jcc = (cc-1)/2;
                                jdd = (dd-1)/2;
                                jee = (ee-1)/2;
                                jff = (ff-1)/2;
                                
                                Zmin = max([jaa+jbb+jee,jee+jcc+jdd,jbb+jcc+jff,jaa+jff+jdd]);
                                Zmax = min([jaa+jbb+jcc+jdd,jaa+jee+jcc+jff,jbb+jee+jdd+jff]);
                                
                                if Zmax>=Zmin
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    
                                    PreNumber = (-1)^(jaa+jbb+jcc+jdd)*sqrt(qNum(ee,q)*qNum(ff,q)...
                                        *(qFactorial(-jaa+jbb+jee,q)*qFactorial(jaa-jbb+jee,q)/qFactorial(jaa+jbb+jee+1,q,jaa+jbb-jee+1))...
                                        *(qFactorial(-jee+jcc+jdd,q)*qFactorial(jee-jcc+jdd,q)/qFactorial(jee+jcc+jdd+1,q,jee+jcc-jdd+1))...
                                        *(qFactorial(-jbb+jcc+jff,q)*qFactorial(jbb-jcc+jff,q)/qFactorial(jbb+jcc+jff+1,q,jbb+jcc-jff+1))...
                                        *(qFactorial(-jaa+jff+jdd,q)*qFactorial(jaa-jff+jdd,q)/qFactorial(jaa+jff+jdd+1,q,jaa+jff-jdd+1)));
                                    
                                    
                                    ValueStep = (-1)^Zmin/qFactorial(Zmin-jaa-jbb-jee,q)/qFactorial(Zmin-jee-jcc-jdd,q)...
                                        /qFactorial(Zmin-jbb-jcc-jff,q)/qFactorial(Zmin-jaa-jff-jdd,q)...
                                        *qFactorial(Zmin+1,q,jaa+jbb+jcc+jdd-Zmin+1)...
                                        /qFactorial(jaa+jee+jcc+jff-Zmin,q)/qFactorial(jbb+jee+jdd+jff-Zmin,q);
                                    
                                    Value = ValueStep;
                                    
                                    for Z = (Zmin+1):Zmax
                                        ValueStep = ValueStep*((-qNum(Z+1,q))...
                                            *qNum(jaa+jbb+jcc+jdd-Z+1,q)*qNum(jaa+jee+jcc+jff-Z+1,q)*qNum(jbb+jee+jdd+jff-Z+1,q)...
                                            /qNum(Z-jaa-jbb-jee,q)/qNum(Z-jee-jcc-jdd,q)/qNum(Z-jbb-jcc-jff,q)/qNum(Z-jaa-jff-jdd,q));
                                        Value = Value+ValueStep;
                                    end
                                    
                                    AccessFValues(Counter,:) = Value*PreNumber;
                                    Counter = Counter+1; 
                                end
                            end
                        end
                    end
                end
                
                AccessFLabelsTemp = AccessFLabels(1:(Counter-1),:);
                [~,~,UniqueFEntries] = unique(AccessFLabelsTemp(:,1:4),'rows');
                
                for jj = 1:max(UniqueFEntries)
                    Loc = find(jj == UniqueFEntries);%this is done in case I need to add extra things
                    MatrixEntries = AccessFValues(Loc,:);
                    SideValues = AccessFLabels(Loc,5:6);%mod this if we need to take multiplicities into account
                    [SideValuesA,~,EntryNumbersA] = unique(SideValues(:,1));
                    [SideValuesB,~,EntryNumbersB] = unique(SideValues(:,2));
                    SparseF = sparse(EntryNumbersA,EntryNumbersB,MatrixEntries,max(max([EntryNumbersA,EntryNumbersB])),max(max([EntryNumbersA,EntryNumbersB])));
                    SparseFInv = SparseF^(-1);
                    [EntryNumbersInvB,EntryNumbersInvA,MatrixEntriesInv] = find(SparseFInv);
                    SideValuesInv = [SideValuesA(EntryNumbersInvA),SideValuesB(EntryNumbersInvB)];
                    
                    [~,UniqueSidesFLocations,UniqueSidesF] = unique([SideValues;SideValuesInv],'rows');
                    
                    %now the important things are the new Values and the
                    %new SideValues
                    
                    %we have a list Loc, which must be extended, first I
                    %assume that all SideValues are unique and will be
                    %selected first so 
                    
                    if sum(UniqueSidesFLocations<=size(SideValues,1))<size(SideValues,1)
                        error('Affirmation Error: UniqueSidesF is wrong, probably because unique is selecting non-SideValues entries as the unique entries')
                    end
                    
                    NewLocations = UniqueSidesFLocations(UniqueSidesFLocations>size(SideValues,1))-size(SideValues,1);
                    
                    if ~isempty(NewLocations)
                        AddedNumber = numel(NewLocations);
                        AccessFLabels(Counter:(Counter+AddedNumber),:) = SideValuesInv(NewLocations,:);
                        
                        Loc = [Loc(:)', Counter-1+(1:AddedNumber)];
                        Counter = Counter+AddedNumber;
                    end
                    
                    %now we need to put in values
                    
                    [~,InverseOriginal] = sort(UniqueSidesF(1:size(SideValues,1)),'ascend');
                    
                    InverseOriginal = [InverseOriginal(:)',(numel(InverseOriginal)+1):numel(UniqueSidesFLocations)];
                    
                    AccessFInvValues(Loc(InverseOriginal(UniqueSidesF((size(SideValues,1)+1):end)))) = MatrixEntriesInv;
                    
                end
                
                
                self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:3
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessRLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRValues = zeros([size(self.Dim,1)^3,1]);
                AccessRInvValues = zeros([size(self.Dim,1)^3,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    if self.Fusion2(cc,aa,bb)>0
                        AccessRLabels(Counter,1:3) = [aa,bb,cc];
                        jaa = (aa-1)/2;
                        jbb = (bb-1)/2;
                        jcc = (cc-1)/2;
                        AccessRValues(Counter,1) = (-1)^(jcc-jaa-jbb)*q^(1/2*(jcc*(jcc+1)-jaa*(jaa+1)-jbb*(jbb+1)));
                        AccessRInvValues(Counter,1) = AccessRValues(Counter,1)';
                        Counter = Counter+1;
                    end
                end
                self.R = {AccessRLabels(1:(Counter-1),:), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                self.FlagBraiding = true;
                
            end
            
            
            if isequal(label,'StMult') %this is for anyonic (non-abelian)
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;3,2,1;3,3,2;4,2,3;5,2,4];
                self.Name = 'StMult';
                
                %fuse with I:
                Entry{1} = eye(5);
                
                %fuse with a:
                Entry{2} = [0,0,1,0,0;...
                            1,3,0,2,0;...
                            0,0,0,0,0;...
                            0,0,2,0,0;...
                            0,0,2,0,3];
                
                %fuse with b:
                Entry{3} = [0,1,0,0,0;...
                            0,0,0,0,0;...
                            1,0,3,0,2;...
                            0,2,0,3,0;...
                            0,2,0,0,0];
                
                %fuse with c:
                Entry{4} = [0,0,0,1,0;...
                            0,2,0,1,0;...
                            0,0,0,0,0;...
                            1,0,3,0,1;...
                            0,0,0,0,1];
                
                %fuse with d:
                Entry{5} = [0,0,0,0,1;...
                            0,0,0,0,0;...
                            0,0,2,0,1;...
                            0,0,0,1,0;...
                            1,3,0,1,0];
                
                
                self.Fusion2 = cat(3,Entry{:});
                
            end
            
            if isequal(label,'SiMult') %this is for anyonic (non-abelian)
                FlagGeneratedSymmetry = true;
                
                self.Dim = [1,1,0;2,2,1;3,2,2;4,1,3];
                self.Name = 'SiMult';
                
                %fuse with I:
                Entry{1} = eye(4);
                
                %fuse with a:
                Entry{2} = [0,0,1,0;...
                            1,2,0,0;...
                            0,0,0,0;...
                            0,0,3,2];
                
                %fuse with b:
                Entry{3} = [0,1,0,0;...
                            0,0,0,0;...
                            1,0,2,0;...
                            0,3,0,2];
                
                %fuse with c:
                Entry{4} = [0,0,0,1;...
                            0,0,0,0;...
                            0,0,0,0;...
                            1,2,2,0];
                
                self.Fusion2 = cat(3,Entry{:});
                
            end
            
            
            
            if length(label)>=3 && isequal(label(1:3),'U1_') %this is for non-finite groups
                FlagGeneratedSymmetry = true;
                
                N = str2num(label(4:end));
                
                if N~=round(N)
                    error('error: We only use integer charges in U1 at this point')
                end
                
                if N <= 0
                    error('error: We need a positive number for U1 at this point')
                end
                
                self.Dim = [(1:(2*N+1))',ones([2*N+1,1]),(-N:N)'];
                self.Name = 'U1';
                self.TrivialIrrep = [N+1];
                self.InverseIrrep = [N:-1:-N]';
                
                for jj = -N:N
                    Fusion2Temp{jj+N+1} = diag(ones([2*N+1-abs(jj),1]),-jj);
                end
                
                self.Fusion2 = cat(3,Fusion2Temp{:});
            end
            
            
            if length(label)>=4 && isequal(label(1:4),'SU2_') %this is for non-finite groups which are non-abelian
                FlagGeneratedSymmetry = true;
                
                if isequal(label(end-[1,0]),'/2')
                    FlagOdd = true;
                else
                    FlagOdd = false;
                end
                N = str2num(label(5:(end-2*FlagOdd)))/2^(FlagOdd);
                
                if ~IsInteger(2*N)
                    error('error: We only use half integer labels for the representations in SU2 at this point')
                end
                
                if N <= 0
                    error('error: We need a positive number for SU2 at this point')
                end
                
                self.Dim = [(1:round(2*N+1))',(1:round(2*N+1))',(0:round(2*N))'/2];
                self.Name = 'SU2';
                self.TrivialIrrep = [1];
                self.InverseIrrep = [1:(2*N+1)]';
                
                for jj = 0:2*N
                    Fusion2Temp{jj+1} = zeros(2*N+1);
                    for kk = 0:2*N
                        Fusion2Temp{jj+1}((abs(jj-kk)+1):2:min((jj+kk+1),round(2*N+1)),kk+1) = ones([floor((min(jj+kk,round(2*N))-abs(jj-kk))/2)+1,1]);
                    end
                end
                self.Fusion2 = cat(3,Fusion2Temp{:});
                
                %now F moves
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:4
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessFLabels = zeros([size(self.Dim,1)^4,6]);
                AccessFValues = zeros([size(self.Dim,1)^4,1]);
                AccessFInvValues = zeros([size(self.Dim,1)^4,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    
                    %check parity
                    if mod(aa+bb+cc+dd,2)==0
                        Listee = (abs(aa-bb)+1):2:((aa+bb)-1);
                        Listee(dd<(abs(Listee-cc)+1)|dd>((Listee+cc)-1))= [];
                        Listff = (abs(cc-bb)+1):2:((cc+bb)-1);
                        Listff(dd<(abs(Listff-aa)+1)|dd>((Listff+aa)-1))= [];
                        
                        if ~isempty(Listee) && ~isempty(Listff)
                            WorkDimInt = [repmat(Listee, [1,numel(Listff)]);reshape(repmat(Listff,[numel(Listee),1]),[1,numel(Listee)*numel(Listff)])];
                            
                            Matrix = zeros([numel(Listee),numel(Listff)]);
                            
                            CounterInt = 1;
                            
                            for WorkInt = WorkDimInt
                                
                                ee = WorkInt(1);
                                ff = WorkInt(2);
                                
                                jaa = (aa-1)/2;
                                jbb = (bb-1)/2;
                                jcc = (cc-1)/2;
                                jdd = (dd-1)/2;
                                jee = (ee-1)/2;
                                jff = (ff-1)/2;
                                
                                
                                Zmin = max([jaa+jbb+jee,jee+jcc+jdd,jbb+jcc+jff,jaa+jff+jdd]);
                                Zmax = min([jaa+jbb+jcc+jdd,jaa+jee+jcc+jff,jbb+jee+jdd+jff]);
                                
                                if Zmax>=Zmin
                                    
                                    PreNumber = (-1)^(jaa+jbb+jcc+jdd)*sqrt(ee*ff...
                                        *(factorial(-jaa+jbb+jee)*factorial(jaa-jbb+jee)*factorial(jaa+jbb-jee)/factorial(jaa+jbb+jee+1))...
                                        *(factorial(-jee+jcc+jdd)*factorial(jee-jcc+jdd)*factorial(jee+jcc-jdd)/factorial(jee+jcc+jdd+1))...
                                        *(factorial(-jbb+jcc+jff)*factorial(jbb-jcc+jff)*factorial(jbb+jcc-jff)/factorial(jbb+jcc+jff+1))...
                                        *(factorial(-jaa+jff+jdd)*factorial(jaa-jff+jdd)*factorial(jaa+jff-jdd)/factorial(jaa+jff+jdd+1)));
                                    
                                    
                                    ValueStep = (-1)^Zmin/factorial(Zmin-jaa-jbb-jee)/factorial(Zmin-jee-jcc-jdd)...
                                        /factorial(Zmin-jbb-jcc-jff)/factorial(Zmin-jaa-jff-jdd)...
                                        *factorial(Zmin+1)/factorial(jaa+jbb+jcc+jdd-Zmin)...
                                        /factorial(jaa+jee+jcc+jff-Zmin)/factorial(jbb+jee+jdd+jff-Zmin);
                                    
                                    Value = ValueStep;
                                    
                                    for Z = (Zmin+1):Zmax
                                        ValueStep = (-1)^Z/factorial(Z-jaa-jbb-jee)/factorial(Z-jee-jcc-jdd)...
                                        /factorial(Z-jbb-jcc-jff)/factorial(Z-jaa-jff-jdd)...
                                        *factorial(Z+1)/factorial(jaa+jbb+jcc+jdd-Z)...
                                        /factorial(jaa+jee+jcc+jff-Z)/factorial(jbb+jee+jdd+jff-Z);
                                        %ValueStep*(-(Z+1)...
                                        %    *(jaa+jbb+jcc+jdd-Z+1)*(jaa+jee+jcc+jff-Z+1)*(jbb+jee+jdd+jff-Z+1)...
                                        %    /(Z-jaa-jbb-jee)/(Z-jee-jcc-jdd)/(Z-jbb-jcc-jff)/(Z-jaa-jff-jdd));
                                        Value = Value+ValueStep;
                                    end
                                    
                                    Matrix(CounterInt) = Value*PreNumber;
                                    CounterInt = CounterInt+1;
                                    
                                end
                            end
                            
                            %now compute inverse and regular
                            InvMatrix = (Matrix^-1).';%need to transpose to get the order correct
                            Keepee = Listee<size(self.Dim,1);
                            Keepff = Listff<size(self.Dim,1);
                            
                            Matrix = Matrix(Keepee,Keepff);
                            InvMatrix = InvMatrix(Keepee,Keepff);
                            
                            ExtraList = WorkDimInt(:,all(WorkDimInt<size(self.Dim,1),1)).';
                            
                            AccessFLabels(Counter-1+(1:size(ExtraList,1)),:) = [repmat([aa,bb,cc,dd],[size(ExtraList,1),1]),ExtraList(:,[2,1])];
                            
                            AccessFValues(Counter-1+(1:size(ExtraList,1)),:) = Matrix(:);
                            AccessFInvValues(Counter-1+(1:size(ExtraList,1)),:) = InvMatrix(:);
                            Counter = Counter+size(ExtraList,1); 
                        end
                    end
                end
                
                self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                
                %{
                AccessFLabels = zeros([size(self.Dim,1)^4,6]);
                AccessFValues = zeros([size(self.Dim,1)^4,1]);
                AccessFInvValues = zeros([size(self.Dim,1)^4,1]);
                Counter = 1;
                for aa = 1:size(self.Dim,1)
                for bb = 1:size(self.Dim,1)
                for cc = 1:size(self.Dim,1)
                for dd = 1:size(self.Dim,1)
                                
                    for ee = 1:size(self.Dim,1)
                        if self.Fusion2(ee,aa,bb)>0 && self.Fusion2(dd,ee,cc)>0
                        for ff = 1:size(self.Dim,1)
                            if self.Fusion2(ff,bb,cc)>0 && self.Fusion2(dd,aa,ff)>0
                                
                                AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                
                                AA = (aa-1)/2;
                                BB = (bb-1)/2;
                                CC = (cc-1)/2;
                                DD = (dd-1)/2;
                                EE = (ee-1)/2;
                                FF = (ff-1)/2;
                                
                                
                                alpha1 = AA+BB+EE;
                                alpha2 = CC+DD+EE;
                                alpha3 = AA+DD+FF;
                                alpha4 = BB+CC+FF;
                                
                                beta1 = AA+BB+CC+DD;
                                beta2 = AA+CC+EE+FF;
                                beta3 = BB+DD+EE+FF;
                                
                                
                                
                                z = max([alpha1,alpha2,alpha3,alpha4]):min([beta1,beta2,beta3]);
                                
                                
                                W = sum(  (-1).^(z+beta1).*factorial(z+1)./factorial(z-alpha1)./factorial(z-alpha2)./factorial(z-alpha3)./factorial(z-alpha4)...
                                    .*(factorial(AA+BB-EE)*factorial(AA-BB+EE)*factorial(-AA+BB+EE)/factorial(AA+BB+EE+1))^(1/2) ...
                                    .*(factorial(CC+DD-EE)*factorial(CC-DD+EE)*factorial(-CC+DD+EE)/factorial(CC+DD+EE+1))^(1/2) ...
                                    .*(factorial(AA+DD-FF)*factorial(AA-DD+FF)*factorial(-AA+DD+FF)/factorial(AA+DD+FF+1))^(1/2) ...
                                    .*(factorial(BB+CC-FF)*factorial(BB-CC+FF)*factorial(-BB+CC+FF)/factorial(BB+CC+FF+1))^(1/2) ...
                                    ./factorial(beta1-z)./factorial(beta2-z)./factorial(beta3-z)  );
                                
                                AccessFValues(Counter,:) = W*(-1)^(AA+BB+CC+DD);
                                AccessFInvValues(Counter,:) = (W*(-1)^(AA+BB+CC+DD))';
                                Counter = Counter+1;
                                
                            end
                        end
                        end
                    end
                end
                end
                end
                end
                self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                %}
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:3
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                AccessRLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRValues = zeros([size(self.Dim,1)^3,1]);
                AccessRInvValues = zeros([size(self.Dim,1)^3,1]);
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    if self.Fusion2(cc,aa,bb)>0
                        
                        AccessRLabels(Counter,:) = [aa,bb,cc];
                        
                        jaa = (aa-1)/2;
                        jbb = (bb-1)/2;
                        jcc = (cc-1)/2;
                        
                        AccessRValues(Counter,:) = (-1)^(jcc-jaa-jbb);
                        AccessRInvValues(Counter,1) = (-1)^(jcc-jaa-jbb);
                        Counter = Counter+1;
                        
                    end
                end
                self.R = {AccessRLabels(1:(Counter-1),:), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                
                self.FlagBraiding = true;
                
                
            end
            
            end
            
            if ~FlagGeneratedSymmetry
                error(['Error: Can''t generate the symmetry - ', label]);
            end
            
            self.UpdateDetails;
            
        end
        
    end
    
    methods(Access = 'private')
        function UpdateDetails(self)
            %now checks for Fusion2 and Dim so that they make sense

            if ~(isinteger(self.Fusion2))
                if ~(isfloat(self.Fusion2))
                    error('Fusion2 should be an integer')
                else
                    Fusion2 = int32(self.Fusion2);
                end
            end
            
            if length(size(self.Fusion2)) == 2
                if ~(size(self.Fusion2) == size(self.Dim,1)*[1,1])
                    error('Fusion2 is not the same size as the number of Dimensions we are considering (and is 2D)')
                end
            elseif ~(size(self.Fusion2) == size(self.Dim,1)*[1,1,1])
                error('Fusion2 is not the same size as the number of Dimensions we are considering')
            end
            
            %now compute the Trivial Irrep from Fusion2
            TempSize = size(self.Dim,1);
            self.TrivialIrrep = [];
            for kk = 1:TempSize
                if (self.Dim(kk,2) == 1)&&( all(all( self.Fusion2(:,:,kk) == eye(TempSize),1),2))
                    self.TrivialIrrep(end+1) = kk;
                end
            end
            
            if isempty(self.TrivialIrrep)
                error('Error: We should have a trivial Irreducible representation')
            end
            if length(self.TrivialIrrep) ~= 1
                error('We should only have one trivial Irreducible representation')
            end
            
            %now compute the InverseIrreps from Fusion2
            self.InverseIrrep = zeros([size(self.Fusion2,1),1]);
            
            for kk = 1:TempSize
                TempInverseIrrep = find(self.Fusion2(self.TrivialIrrep, :, kk));
                if numel(TempInverseIrrep)~=1
                    if numel(TempInverseIrrep) == 0
                        warning(['This symmetry ', self.Name, ' has no inverse Irreps for weight ' num2str(self.Dim(kk,2:end))])
                    else
                        error(['This symmetry ', self.Name, ' has multiple inverse Irreps for weight ' num2str(self.Dim(kk,2:end))])
                    end
                end
                self.InverseIrrep(kk) = TempInverseIrrep;
            end
            
            %check if we are self dual
            if self.InverseIrrep == (1:length(self.InverseIrrep))'
                self.FlagSelfDual = true;
            end
            
            %check if we are non-abelian
            if any(any(sum(self.Fusion2,1)>1,2),3)||any(self.Dim(:,2)>1)
                self.FlagNonAbelian = true;
            end
            
            %check if we have any multiplicities
            self.MaxCharges = size(self.Dim,1);
            self.MaxMultiplicities = max(max(max(self.Fusion2)));
            if self.MaxMultiplicities>1
                self.FlagMultiplicities = true;
            end
            
            
            if (self.FlagBraiding) || (self.FlagNonAbelian)||true
                
                
                if ~self.FlagBraiding
                    AccessRLabels = zeros([size(self.Dim,1)^3*self.MaxMultiplicities,5]);
                    AccessRValues = zeros([size(self.Dim,1)^3*self.MaxMultiplicities,1]);
                    AccessRInvValues = zeros([size(self.Dim,1)^3*self.MaxMultiplicities,1]);
                    
                    
                    WorkDim = 1:size(self.Dim,1);
                    for kk = 2:3
                        WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                    end
                    
                    Counter = 1;
                    
                    for Work = WorkDim
                        aa = Work(1);
                        bb = Work(2);
                        cc = Work(3);
                        for mm = 1:self.Fusion2(cc,aa,bb)
                            AccessRLabels(Counter,:) = [aa,bb,cc,mm,mm];
                            AccessRValues(Counter,:) = 1;
                            AccessRInvValues(Counter,1) = 1;
                            Counter = Counter+1;
                        end
                    end
                    
                    if self.FlagMultiplicities
                        self.R = {AccessRLabels(1:(Counter-1),:), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                    else
                        self.R = {AccessRLabels(1:(Counter-1),1:3), AccessRValues(1:(Counter-1),:), AccessRInvValues(1:(Counter-1),:)};
                    end
                end
                
                if ~self.FlagNonAbelian
                    AccessFLabels = zeros([size(self.Dim,1)^6,6]);
                    AccessFValues = zeros([size(self.Dim,1)^6,1]);
                    AccessFInvValues = zeros([size(self.Dim,1)^6,1]);
                    
                    %Note that in this case MaxMultiplicities = 1
                    
                    if self.MaxMultiplicities>1
                        error('Affimation Error: We should be multiplicity free here');
                    end
                    
                    Counter = 1; 
                    WorkDim = 1:size(self.Dim,1);
                    for kk = 2:4
                        WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                    end
                    
                    for Work = WorkDim
                        aa = Work(1);
                        bb = Work(2);
                        cc = Work(3);
                        dd = Work(4);
                        
                        Possible_ee = find(self.Fusion2(:,aa,bb));
                        Possible_ff = find(self.Fusion2(:,bb,cc));
                        if ~isempty(Possible_ee)&&~isempty(Possible_ff)
                            PossibleDim = [repmat(Possible_ee(:),[numel(Possible_ff),1]),reshape(repmat(Possible_ff(:)',[numel(Possible_ee),1]),[numel(Possible_ff)*numel(Possible_ee),1])];
                            for Work2 = PossibleDim'
                                
                                ee = Work2(1);
                                ff = Work2(2);
                                if self.Fusion2(ee,aa,bb)==1 && self.Fusion2(dd,ee,cc)==1 &&...
                                            self.Fusion2(ff,bb,cc)==1 && self.Fusion2(dd,aa,ff)==1
                                    AccessFLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                    AccessFValues(Counter,:) = 1;
                                    AccessFInvValues(Counter,1) = 1;
                                    Counter = Counter+1;
                                end
                            end
                        end
                    end
                    
                    self.F = {AccessFLabels(1:(Counter-1),:), AccessFValues(1:(Counter-1),:), AccessFInvValues(1:(Counter-1),:)};
                end
                
                %now make the FS indicator
                LInverseIrrep = self.InverseIrrep;
                LTrivialIrrep = self.TrivialIrrep;
                FLabels = self.F{1};
                FMoves = self.F{2};
                FInvMoves = self.F{3};
                RLabels = self.R{1};
                RMoves = self.R{2};
                RInvMoves = self.R{3};
                for kk = 1:size(self.Dim,1)
                    Location = find(all(FLabels==repmat([kk,LInverseIrrep(kk),kk,kk,LTrivialIrrep,LTrivialIrrep],[size(FLabels,1),1]),2));
                    if numel(Location) ~=1
                        error('Affirmation Error: We should have a single location when computing Forbenius Schur indicator')
                    end
                    self.Frobenius_Schur_Indicator(kk) = FMoves(Location)*self.Dim(kk,2);
                end
                
                
                for kk = 1:size(self.Dim,1)
                    Location = find(all(RLabels==repmat([kk,LInverseIrrep(kk),LTrivialIrrep],[size(RLabels,1),1]),2));
                    if numel(Location) ~=1
                        error('Affirmation Error: We should have a single location when computing Forbenius Schur indicator')
                    end
                    self.Round(kk) = RInvMoves(Location);
                    %self.RoundInv(kk) = RMoves(Location);
                end
                
                
                %now make the RCCW (A) and RCW (B):
                AccessRCWLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRCWValues = zeros([size(self.Dim,1)^3,1]);
                AccessRCWInvValues = zeros([size(self.Dim,1)^3,1]);
                
                AccessRCCWLabels = zeros([size(self.Dim,1)^3,3]);
                AccessRCCWValues = zeros([size(self.Dim,1)^3,1]);
                AccessRCCWInvValues = zeros([size(self.Dim,1)^3,1]);
                
                Counter = 1;
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:3
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    if self.Fusion2(cc,aa,bb)>0
                    
                        AccessRCCWLabels(Counter,:) = [aa,bb,cc];
                        AccessRCWLabels(Counter,:) = [aa,bb,cc];
                        
                        Location = find(all(FLabels == repmat([LInverseIrrep(aa),aa,bb,bb,cc,LTrivialIrrep],[size(FLabels,1),1]),2));
                        if numel(Location)>1;
                            error('Error: can''t find location for F for RCCW')
                        end
                        
                        if ~isempty(Location)
                            AccessRCCWValues(Counter,:) = sqrt(self.Dim(aa,2)*self.Dim(bb,2)/self.Dim(cc,2)) * self.Frobenius_Schur_Indicator(aa)'*FMoves(Location)';
                            AccessRCCWInvValues(Counter,:) = (sqrt(self.Dim(aa,2)*self.Dim(bb,2)/self.Dim(cc,2)) * self.Frobenius_Schur_Indicator(aa)'*FMoves(Location)')^-1;
                        else
                            AccessRCCWValues(Counter,:) = 0;
                            AccessRCCWInvValues(Counter,:) = 0;
                        end
                        
                        
                        
                        Location = find(all(FLabels == repmat([aa,bb,LInverseIrrep(bb),aa,LTrivialIrrep,cc],[size(FLabels,1),1]),2));
                        if numel(Location)>1;
                            error('Error: can''t find location for F for RCCW')
                        end
                        
                        if ~isempty(Location)
                            AccessRCWValues(Counter,:) = sqrt(self.Dim(aa,2)*self.Dim(bb,2)/self.Dim(cc,2))*FMoves(Location);
                            AccessRCWInvValues(Counter,:) = (sqrt(self.Dim(aa,2)*self.Dim(bb,2)/self.Dim(cc,2))*FMoves(Location))^-1;
                        else
                            AccessRCWValues(Counter,:) = 0;
                            AccessRCWInvValues(Counter,:) = 0;
                        end
                        
                        %note: the dimension out the front corrects
                        %for the fact that FMoves'*FMoves<1 by def
                        %so we expect it to be the inverse of the
                        %dimensions (note that cc is always smaller
                        %then aa * bb by def (don't understand
                        %where the cc comes from).
                        Counter = Counter+1;
                    end
                end
                self.RCW = {AccessRCWLabels(1:(Counter-1),:), AccessRCWValues(1:(Counter-1),:), AccessRCWInvValues(1:(Counter-1),:)};
                self.RCCW = {AccessRCCWLabels(1:(Counter-1),:), AccessRCCWValues(1:(Counter-1),:), AccessRCCWInvValues(1:(Counter-1),:)};
                
                %now make the B (Braiding) and FSq:
                AccessBLabels = zeros([size(self.Dim,1)^3,6]);
                AccessBValues = zeros([size(self.Dim,1)^3,1]);
                AccessBInvValues = zeros([size(self.Dim,1)^3,1]);
                
                AccessFSqLabels = zeros([size(self.Dim,1)^3,6]);
                AccessFSqValues = zeros([size(self.Dim,1)^3,1]);
                AccessFSqInvValues = zeros([size(self.Dim,1)^3,1]);
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:5
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    ee = Work(5);
                    
                    if self.Fusion2(dd,ee,bb)>0 && self.Fusion2(ee,aa,cc)>0
                        for ff = 1:size(self.Dim,1)
                            if self.Fusion2(ff,aa,bb)>0 && self.Fusion2(dd,ff,cc)>0
                                AccessBValues(Counter,:) = 0;
                                AccessBInvValues(Counter,:) = 0;
                                AccessBLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                
                                for gg = 1:size(self.Dim,1)
                                    LocationF1 = find(all(FLabels == repmat([aa,cc,bb,dd,gg,ee],[size(FLabels,1),1]),2));
                                    if numel(LocationF1)>1;
                                        error('Error: can''t find location for F for RCCW')
                                    end
                                    LocationF2 = find(all(FLabels == repmat([aa,bb,cc,dd,gg,ff],[size(FLabels,1),1]),2)); %problem???
                                    if numel(LocationF2)>1;
                                        error('Error: can''t find location for F for RCCW')
                                    end
                                    LocationR1 = find(all(RLabels == repmat([bb,cc,gg],[size(RLabels,1),1]),2));
                                    if numel(LocationR1)>1;
                                        error('Error: can''t find location for F for RCCW')
                                    end
                                    
                                    if ~isempty(LocationF1)&&~isempty(LocationF2)&&~isempty(LocationR1)
                                        AccessBValues(Counter,:) = AccessBValues(Counter,:) + FMoves(LocationF1)*RMoves(LocationR1)*FInvMoves(LocationF2)'; %problem???
                                        AccessBInvValues(Counter,:) = AccessBInvValues(Counter,:) + FInvMoves(LocationF1)*RInvMoves(LocationR1)*FMoves(LocationF2)';
                                    end
                                end
                                Counter = Counter+1;
                            end
                        end
                    end
                end
                self.B = {AccessBLabels(1:(Counter-1),:), AccessBValues(1:(Counter-1),:), AccessBInvValues(1:(Counter-1),:)};
                
                
                
                WorkDim = 1:size(self.Dim,1);
                for kk = 2:5
                    WorkDim = [repmat(WorkDim,[1,size(self.Dim,1)]);reshape(repmat(1:size(self.Dim,1),[size(WorkDim,2),1]),[1,size(self.Dim,1)*size(WorkDim,2)])];
                end
                
                Counter = 1;
                for Work = WorkDim
                    aa = Work(1);
                    bb = Work(2);
                    cc = Work(3);
                    dd = Work(4);
                    ee = Work(5);
                    if self.Fusion2(aa,cc,ee)>0 && self.Fusion2(dd,ee,bb)>0
                        for ff = 1:size(self.Dim,1)
                            if self.Fusion2(ff,aa,bb)>0 && self.Fusion2(ff,cc,dd)>0
                                
                                AccessFSqLabels(Counter,:) = [aa,bb,cc,dd,ff,ee];
                                
                                LocationF1 = find(all(FLabels == repmat([LInverseIrrep(cc),aa,bb,dd,ff,ee],[size(FLabels,1),1]),2));
                                if numel(LocationF1)>1;
                                    error('Error: can''t find location for F for RCCW')
                                end
                                LocationA1 = find(all(AccessRCCWLabels == repmat([LInverseIrrep(cc),aa,ee],[size(AccessRCCWLabels,1),1]),2));
                                if numel(LocationA1)>1;
                                    error('Error: can''t find location for F for RCCW')
                                end
                                LocationA2 = find(all(AccessRCCWLabels == repmat([LInverseIrrep(cc),ff,dd],[size(AccessRCCWLabels,1),1]),2));
                                if numel(LocationA2)>1;
                                    error('Error: can''t find location for F for RCCW')
                                end
                                
                                if ~isempty(LocationF1)&&~isempty(LocationF2)&&~isempty(LocationR1)
                                    AccessFSqValues(Counter,:) = FMoves(LocationF1)*AccessRCCWInvValues(LocationA1)*AccessRCCWValues(LocationA2);
                                    AccessFSqInvValues(Counter,:) = FInvMoves(LocationF1)*AccessRCCWValues(LocationA1)*AccessRCCWInvValues(LocationA2);
                                else
                                    AccessFSqValues(Counter,:) = 0;
                                    AccessFSqInvValues(Counter,:) = 0;
                                end
                                Counter = Counter+1;
                            end
                        end
                    end
                end
                self.FSq = {AccessFSqLabels(1:(Counter-1),:), AccessFSqValues(1:(Counter-1),:), AccessFSqInvValues(1:(Counter-1),:)};
                
            end 
            
        end
    end
    
end

function bool = IsInteger(A)
    bool = isnumeric(A);
    bool = bool&~any(abs(A(:)-round(A(:)))>10^-14);
end

function Out = qNum(N,q)
    Out = (q^(N/2)-q^(-N/2))/(q^(1/2)-q^(-1/2));
end

function Out = qFactorial(N,q,Low)
    if nargin<3
        Low = 1;
    end
    
    if N == 0
        Out = 1;
        return;
    end
    
    if Low<1
        error('Error: Low should be greater then or equal to 1')
    end
    
    Value = 1;
    for nn = Low:N
        Value = Value*qNum(nn,q);
    end
    Out = Value;
    
end

function [C,IA,IC] = winunique(A, winSz)
    if nargin<2
        winSz = 10^2;
    end
    
    if size(A,1)<winSz
        [C,IA,IC] = unique(A,'rows');
    else 
        nRows = size(A,1);
        C = zeros(0,size(A,2));
        for k=1:winSz:nRows
            [C2, IA2, IC2] = unique([C; A(k:min(k+winSz-1,end),:)], 'rows');
            
            IA2 = IA2+k-1-size(C,1);
            if numel(C)>0
                IA2(IA2<=(size(C,1)+k-1-size(C,1))) = IA(IA2(IA2<=(size(C,1)+k-1-size(C,1)))-k+1+size(C,1));
                IA = IA2;
                
                IC = IC2(IC);
                IC = [IC;IC2((size(C,1)+1):end)];
            else
                IC = IC2;
                IA = IA2;
            end
            
            C = C2;
        end
        %[~, IA] = unique(I);
        %C = A(IA, :);
    end
end

%Based on code taken from:
%https://stackoverflow.com/questions/29829087/octave-matlab-determine-unique-rows-in-very-large-matrix


%{
function [LeftIndiciesOut, RightIndiciesOut, EntriesOut] = HackSparse(LeftIndiciesInA, RightIndiciesInA, EntriesInA, LeftIndiciesInB, RightIndiciesInB, EntriesInB)
    
    ZeroBound = 10^-14;
    
    %don't need first locaiton of a unique value
    [UniqueLeft,~,UniqueValueNumberLeft] = unique(LeftIndiciesInA);
    [UniqueRight,~,UniqueValueNumberRight] = unique(RightIndiciesInB);
    
    [~,~,UniqueValueNumberMid] = unique([RightIndiciesInA;LeftIndiciesInB]);
    
    UniqueValueNumberMidA = UniqueValueNumberMid(1:numel(RightIndiciesInA));
    UniqueValueNumberMidB = UniqueValueNumberMid((numel(RightIndiciesInA)+1):end);
    
    ForLoopList = [repmat(reshape(1:numel(UniqueLeft),[numel(UniqueLeft),1]),[1,1,numel(UniqueRight)]),...
        repmat(reshape(1:numel(UniqueRight),[1,1,numel(UniqueRight)]),[numel(UniqueLeft),1,1])];
    
    ForLoopList = reshape(permute(ForLoopList,[1,3,2]),[numel(UniqueLeft)*numel(UniqueRight),2])';
    
    LeftIndiciesOut = zeros([0,1]);
    RightIndiciesOut = zeros([0,1]);
    EntriesOut = zeros([0,1]);
    
    for ff = ForLoopList
        
        EntriesPartsA = EntriesInA(UniqueValueNumberLeft == ff(1));
        MidIndiciesA = UniqueValueNumberMidA(UniqueValueNumberLeft == ff(1));
        
        EntriesPartsB = EntriesInB(UniqueValueNumberRight == ff(2));
        MidIndiciesB = UniqueValueNumberMidB(UniqueValueNumberRight == ff(2));
        
        TotalSum = 0;
        for kk = 1:numel(MidIndiciesA)
            Temp = MidIndiciesB==MidIndiciesA(kk);
            
            if sum(Temp)>1
                error('Affirmation Error: there should not be more then one pair in the entry')
            elseif sum(Temp)>0
                TotalSum = TotalSum + EntriesPartsA(kk)*EntriesPartsB(Temp);
            end
        end
        
        
        if abs(TotalSum) > ZeroBound
            LeftIndiciesOut = [LeftIndiciesOut;UniqueLeft(ff(1))];
            RightIndiciesOut = [RightIndiciesOut;UniqueRight(ff(2))];
            EntriesOut = [EntriesOut;TotalSum];
        end
        
    end
    
    
    
end
%}