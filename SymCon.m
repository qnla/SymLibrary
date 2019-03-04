function TensorsOut = SymCon(Tensors, Legs, ContractionOrder, LegOrder, FlagNumber, Crossings)
%This function is suppose to act as an abstraction function between the
%Acore function and the MERA function, this wrapping means that I can also
%swap easily between the ACore and my functions by swapping the code
%around, the 
%
%

MaxChoices = 10000;

FlagAllowBraidingOldTrue = 1==0;
FlagPerformPlanarStuff = 1==1;
FlagSpeedUpAbelian = true;
 
MinNonZero = 10^-14;

if nargin<4
    LegOrder = {[]};
end

if nargin<5
    FlagNumber = 0;
else
    if ~isnumeric(FlagNumber)
        FlagNumber = 0;
    elseif ~IsInteger(FlagNumber)
        FlagNumber = 0;
    elseif FlagNumber<0
        FlagNumber = 0;
    end
end

if nargin<6
    Crossings = [];
end

%first check that the tensors are Acore, and record the dimensions

if ~iscell(Tensors)
    error('Tensors should be a 1xN cell of SymTensor valued entries')
end
if numel(Tensors)<1
    error('There should be some Tensor in the list of Tensors')
end
Tensors = Tensors(:);



if isempty(Tensors{1})
    Tensors(1) = [];
    FlagSingleNotCell = true;
    if numel(Tensors)<1
       error('There should be some Tensor in the list of Tensors after the empty entry')
    end
else
    FlagSingleNotCell = false;
end

NumberLegsTensor  = zeros([1,numel(Tensors)]);

ExternalChargeTensorDetails = cell(size(Tensors));
InternalChargeTensorDetails = cell(size(Tensors));
MultiplicitiesInternalTensorDetails = cell(size(Tensors));
ChargeLegDimensionsTensorDetails = cell(size(Tensors));
ChargeSideTensorDetails = cell(size(Tensors));
ChargeSideIntTensorDetails = cell(size(Tensors));
StructureTensorDetails = cell(size(Tensors));
ChargeDirectionsTensorDetails = cell(size(Tensors));
CyclicLegsBraidedTensorDetails = cell(size(Tensors));



for kk = 1:length(Tensors)
    if ~isa(Tensors{kk},'SymTensor')
        error('Entry %i in Tensors is not an SymTensor valued object', kk)
    end
    NumberLegsTensor(kk) = nLegs(Tensors{kk});
    
    if kk == 1
        SymmetryUsed = Tensors{kk}.Symmetry;
    else
        if ~isequal(SymmetryUsed,Tensors{kk}.Symmetry)
            error('All Tensors should have the same symmetry (in fact the same instance of the symmetry)')
        end
    end
    Tensors{kk}.SortLabels;
    
    ExternalChargeTensorDetails{kk} = Tensors{kk}.ChargeLabelsExternal;
    InternalChargeTensorDetails{kk} = Tensors{kk}.ChargeLabelsInternal;
    MultiplicitiesInternalTensorDetails{kk} = Tensors{kk}.MultiplicitiesInternal;
    ChargeLegDimensionsTensorDetails{kk} = Tensors{kk}.ChargeLegDimensions;
    ChargeSideTensorDetails{kk} = Tensors{kk}.ChargeSide;
    ChargeSideIntTensorDetails{kk} = Tensors{kk}.ChargeSideInt;
    StructureTensorDetails{kk} = Tensors{kk}.Structure;
    ChargeDirectionsTensorDetails{kk} = Tensors{kk}.ChargeDirections;
    CyclicLegsBraidedTensorDetails{kk} = Tensors{kk}.CyclicLegsBraided;
end





if FlagNumber>0
    FlagMustUsePlanar = SymmetryUsed.IsNonAbelian || SymmetryUsed.IsBraiding;
    
    FlagLoadSaveBroad = SymmetryUsed.CheckValidLoad(FlagNumber,  Legs, ContractionOrder, LegOrder,ChargeSideTensorDetails,StructureTensorDetails,ChargeDirectionsTensorDetails);
    FlagLoadSavePlanar = FlagLoadSaveBroad&&(FlagMustUsePlanar);
    FlagNewBroadDetails = ~FlagLoadSaveBroad;
    FlagNewPlanarDetails = FlagNewBroadDetails&&(FlagPerformPlanarStuff||FlagMustUsePlanar);
    
    FlagLoadSaveMid = SymmetryUsed.CheckValidLoadMid(FlagNumber,ExternalChargeTensorDetails,InternalChargeTensorDetails,MultiplicitiesInternalTensorDetails) ...
        && ~isequal(SymmetryUsed.getSymName,'none');
    FlagLoadSavePlanarMid = FlagLoadSaveMid && (FlagMustUsePlanar);
    FlagNewMidDetails = ~FlagLoadSaveMid  &&  ~isequal(SymmetryUsed.getSymName,'none');
    
    FlagLoadSaveFine = SymmetryUsed.CheckValidLoadFine(FlagNumber, ChargeLegDimensionsTensorDetails);
    FlagNewFineDetails = ~FlagLoadSaveFine;
    
    
    FlagSaveAtEnd = ~FlagLoadSaveBroad;
    FlagSavePlanarAtEnd = ~FlagLoadSaveBroad&&(FlagMustUsePlanar);
    FlagSaveMidAtEnd = ((~FlagSaveAtEnd&&FlagNewMidDetails)||~FlagLoadSaveMid)&&~isequal(SymmetryUsed.getSymName,'none');
    FlagSaveMidPlanarAtEnd = (FlagSaveMidAtEnd || FlagSaveAtEnd)&&FlagMustUsePlanar;
    FlagSaveFineAtEnd = ~FlagSaveMidAtEnd||~FlagLoadSaveFine;
    
    if FlagSaveAtEnd
        SavedCheckDetails = {Legs, ContractionOrder, LegOrder, ChargeSideTensorDetails,StructureTensorDetails,ChargeDirectionsTensorDetails};
        SavedCheckMidDetails = {ExternalChargeTensorDetails,InternalChargeTensorDetails,MultiplicitiesInternalTensorDetails};
        SavedCheckFineDetails = {ChargeLegDimensionsTensorDetails};
    end
    
    if FlagSaveMidAtEnd
        SavedCheckMidDetails = {ExternalChargeTensorDetails,InternalChargeTensorDetails,MultiplicitiesInternalTensorDetails};
        SavedCheckFineDetails = {ChargeLegDimensionsTensorDetails};
    end
    
    if FlagSaveFineAtEnd
        SavedCheckFineDetails = {ChargeLegDimensionsTensorDetails};
    end
else
    FlagMustUsePlanar = SymmetryUsed.IsNonAbelian || SymmetryUsed.IsBraiding;
    
    FlagLoadSaveBroad = false;
    FlagLoadSavePlanar = false;
    FlagNewBroadDetails = true;
    FlagNewPlanarDetails = (FlagPerformPlanarStuff||FlagMustUsePlanar);
    
    FlagLoadSaveMid = false;
    FlagLoadSavePlanarMid = false;
    FlagNewMidDetails = ~isequal(SymmetryUsed.getSymName,'none');
    
    FlagLoadSaveFine = false;
    FlagNewFineDetails = true;
    
    
    FlagSaveAtEnd = false;
    FlagSavePlanarAtEnd = false;
    FlagSaveMidAtEnd = false;
    FlagSaveMidPlanarAtEnd = false;
    FlagSaveFineAtEnd = false;
end

%disp(['broad: ', num2str(FlagLoadSaveBroad),' & Mid: ', num2str(FlagLoadSaveMid)])

if FlagNewBroadDetails
                    
    [LegOrder,FlagUserDidNotChooseLegOrder,UniqueLegsNeg, TotalTensors,TensorCopies] =...
    Broad_1_ComputingWorkingLegs(LegOrder,Legs,ChargeSideTensorDetails, ChargeDirectionsTensorDetails,StructureTensorDetails,NumberLegsTensor);
    
    
    [UniqueLegsPos, ModifiedListIncluded,TracingTreeContractions, TracingTreeIncluded, TracingTreeTriads,...
    MaxNumberEdges, AtomicTensorTracing, AtomicTensorOutputTracing, TracedLegsContractions, NegativesTopTensors, NegativeReplacements,...
    FullContractionClosedFinalTriad, FullContractionClosedFinalContractedLegs] = ...
    Broad_2_ComputingTreeStructure(TensorCopies, ContractionOrder, TotalTensors, Legs, Crossings);
    
    
    [NeedFullContractionClosed, EnvironmentLabelsForFullContractions, RequiredAtomicTensors, EnvironmentLabelsForAtomicTensors, NeedCenterTriadForAbove,...
    EnvironmentLabelsForCenterTriadAbove, NeedTriadsForAbove, EnvironmentLabelsForTriadsAbove, NeedTriadsForBelow, EnvironmentLabelsForTriadsBelow] = ...
    Broad_3_ComputingEnvironments(LegOrder, TotalTensors, NegativesTopTensors, TracingTreeTriads, TracingTreeIncluded, ModifiedListIncluded );
    
    
    [OldTrue, OldTrueLocations, ActionTotalTensors, ActionTensorsLocations, ActualLegsTotal, DropTensor, ActionTypeTotal, LegsContractedTotal,EnvironmentLegs] = ...
    Broad_4_ComputingContractionOrder(NeedFullContractionClosed, EnvironmentLabelsForFullContractions, RequiredAtomicTensors, EnvironmentLabelsForAtomicTensors, ...
    NeedCenterTriadForAbove, EnvironmentLabelsForCenterTriadAbove, NeedTriadsForAbove, EnvironmentLabelsForTriadsAbove, NeedTriadsForBelow, EnvironmentLabelsForTriadsBelow, ...
    LegOrder, TotalTensors, TracingTreeTriads, TracingTreeContractions, MaxNumberEdges, ModifiedListIncluded, AtomicTensorTracing, AtomicTensorOutputTracing, TracedLegsContractions, ...
    Legs, NegativeReplacements,FullContractionClosedFinalTriad, FullContractionClosedFinalContractedLegs);
    
    
    [OldTrueChargeSidesInt, OldTrueBraiding, OldTrueBraidingDirection, OldTrueStructure, OldTrueChargeDirections,...
        ActionTensorsPermute,ActionNumberContractedLegs,OldTrueChargeSides,OldTruePermute,ActionTensorsStructure,ActionTensorsChargeDirections] = ...
    Broad_5_ComputingOldTrue(LegOrder, EnvironmentLegs, StructureTensorDetails, ChargeDirectionsTensorDetails, ActualLegsTotal, LegsContractedTotal, ...
    ActionTotalTensors, ActionTypeTotal, ActionTensorsLocations, OldTrue, OldTrueLocations, TotalTensors);
    
    
end

if FlagNewPlanarDetails
    
    if ~isempty(UniqueLegsNeg)
        UniqueLegsPosOld = UniqueLegsPos;
        
        RenameNeg = [UniqueLegsNeg;max(UniqueLegsPos)+(1:numel(UniqueLegsNeg))];
        UniqueLegsPos = [UniqueLegsPos, max(UniqueLegsPos)+(1:numel(UniqueLegsNeg))];
        
        KeepNegs = [];
        NegChargeSides=zeros([1,max(UniqueLegsNeg)]);
        for kk =1:numel(ChargeSideTensorDetails)
            NegLocs = (Legs{kk} <0);
            
            if ~isempty(NegLocs)
                FoundNegs = -Legs{kk}(NegLocs);
                KeepNegs = [KeepNegs,-FoundNegs];
                NegChargeSides(FoundNegs) = ChargeSideTensorDetails{kk}(NegLocs);
            end
            
        end
        
        KeepNegs = sort(KeepNegs,'descend');
        
        if ~isequal(KeepNegs,UniqueLegsNeg)
            error('Affirmation Error: KeepNegs should be the same as UniqueLegsNeg')
        end
        
        
        %ERROR HERA, assumes they are ordered
        NegChargeSides = NegChargeSides(-KeepNegs);
        CyclicLegsNeg = RenameNeg(2,NegChargeSides>0);
        CyclicLegsNeg = [RenameNeg(2,NegChargeSides<0),CyclicLegsNeg(end:-1:1)];
        
    end
    
    
    %first work out the loops
    
    if ~isempty(UniqueLegsNeg)
        CyclicLegsBraidedLabelsUser = cell([TotalTensors(end)+1,1]);
        CyclicLegsBraidedUsed = cell([TotalTensors(end)+1,1]);
        CyclicLegsBraidedUsedNew = cell([TotalTensors(end),1]);
        CyclicLegsBraidedTopBottomSide = cell([TotalTensors(end)+1,1]);
        
        
        for kk = 1:numel(Legs)
            LegsMod{kk} = Legs{kk};
            for ll = 1:size(RenameNeg,2)
                
                LegsMod{kk}(LegsMod{kk}==RenameNeg(1,ll)) = RenameNeg(2,ll);
                
            end
        end
        NumberTensors = numel(CyclicLegsBraidedTensorDetails);
        CyclicLegsBraidedTensorDetailsOld = CyclicLegsBraidedTensorDetails;
        
    else
        CyclicLegsBraidedLabelsUser = cell([TotalTensors(end),1]);
        CyclicLegsBraidedUsed = cell([TotalTensors(end),1]);
        CyclicLegsBraidedUsedNew = cell([TotalTensors(end),1]);
        CyclicLegsBraidedTopBottomSide = cell([TotalTensors(end),1]);
        NumberTensors = numel(CyclicLegsBraidedTensorDetails);
        LegsMod = Legs;
    end
    
    for kk = 1:(numel(TotalTensors)-1)
        NumberLegs((TotalTensors(kk)+1):TotalTensors(kk+1)) = size(LegsMod{kk},2);
        CyclicLegsBraidedLabelsUser((TotalTensors(kk)+1):TotalTensors(kk+1)) = mat2cell(LegsMod{kk}(:,CyclicLegsBraidedTensorDetails{kk}),ones([size(LegsMod{kk},1),1]));
        CyclicLegsBraidedUsed((TotalTensors(kk)+1):TotalTensors(kk+1)) = mat2cell(false(size(LegsMod{kk})),ones([size(LegsMod{kk},1),1]));
        %if isempty(UniqueLegsNeg)||kk ~=numel(CyclicLegsBraidedTensorDetails)
            CyclicLegsBraidedUsedNew((TotalTensors(kk)+1):TotalTensors(kk+1)) = mat2cell(Legs{kk}(:,CyclicLegsBraidedTensorDetails{kk})<0,ones([size(Legs{kk},1),1]));
        %end
    end
    
    for kk = 1:NumberTensors
        CyclicLegsBraidedTopBottomSide{kk} = ChargeSideTensorDetails{kk}(CyclicLegsBraidedTensorDetails{kk});
    end
    ExternalTensorNumber = 0;
    if ~isempty(UniqueLegsNeg)
        CyclicLegsBraidedLabelsUser{end} = CyclicLegsNeg;
        CyclicLegsBraidedTopBottomSide{end} = -NegChargeSides;
        NumberLegs(TotalTensors(end)+1) = numel(CyclicLegsNeg);
        CyclicLegsBraidedUsed{end} = false(size(CyclicLegsNeg));
        ExternalTensorNumber = numel(CyclicLegsBraidedUsed);
    end
    
    CyclicLegsBraidedFirstSide = CyclicLegsBraidedUsed;
    CyclicLegsBraidedExternalLeg = CyclicLegsBraidedUsed; %we will want to store
    CyclicLegsBraidedFirstSideNew = CyclicLegsBraidedUsedNew;
    CyclicLegsBraidedExternalLegNew = CyclicLegsBraidedUsedNew; %we will want to store
    
    RotateTensor = zeros([numel(CyclicLegsBraidedTensorDetails),1]);
    
    %convert to numbers:
    if ~isempty(UniqueLegsNeg)
        LegsUserLabels = cell2mat([CyclicLegsBraidedLabelsUser(:)',RenameNeg(2,:)]);
        [Original,~,LegsCompLabels] = unique(LegsUserLabels(LegsUserLabels>0));
        RenameNeg(2,:) = LegsCompLabels(end+1-(size(RenameNeg,2):-1:1));
        LegsCompLabels(end+1-(1:size(RenameNeg,2))) = [];
        LegsUserLabels(end+1-(1:size(RenameNeg,2))) = [];
    else
        RenameNeg = zeros([2,0]);
        LegsUserLabels = cell2mat(CyclicLegsBraidedLabelsUser(:)');
        [Original,~,LegsCompLabels] = unique(LegsUserLabels(LegsUserLabels>0));
    end
    
    LegsTranslateOldNew = [reshape(Original,[numel(Original),1]),(1:numel(Original))'];
    %leave this in for debugging
    
    LegsUserLabels(LegsUserLabels>0) = LegsCompLabels;
    CyclicLegsBraidedLabelsUser = mat2cell(LegsUserLabels, 1, NumberLegs);
    CyclicLegsBraidedLabelsUserOld = CyclicLegsBraidedLabelsUser;
    
    NumberLegsContracted = zeros([1,numel(LegsContractedTotal)]);
    for kk = 1:numel(LegsContractedTotal)
        NumberLegsContracted(kk) = numel(LegsContractedTotal{kk});
    end
    
    LegsContractedUserLabels = [Original,cell2mat(LegsContractedTotal(:)')];
    [~,~,LegsContractedUserLabels] = unique(LegsContractedUserLabels);
    LegsContractedTotal = mat2cell(reshape(LegsContractedUserLabels((numel(Original)+1):end),[1,sum(NumberLegsContracted)]), 1, NumberLegsContracted);
    
    if any(LegsContractedUserLabels > numel(Original))
        error('Affirmation Error: a label appeared which didn''t appear on CyclicLegsBraidedLabelsUser')
    end
    
    %now we can assume that they are labelled 1 to N
    
    LegsCounted = zeros([numel(UniqueLegsPos),1]);
    LegLocations = zeros([numel(UniqueLegsPos),2]);
    LegSides = LegLocations;
    
    LegLocationsNeg = zeros([numel(UniqueLegsNeg),1]);
    LegSidesNeg = LegLocationsNeg;
    
    for kk = 1:numel(CyclicLegsBraidedLabelsUser)
        CurrentSides = CyclicLegsBraidedTopBottomSide{kk};
        CurrentLegs = CyclicLegsBraidedLabelsUser{kk};
        for ll = 1:numel(CurrentLegs)
            if CurrentLegs(ll)>0
                if LegLocations(CurrentLegs(ll),1)~=0
                    LegLocations(CurrentLegs(ll),2) = kk;
                    LegSides(CurrentLegs(ll),2) = CurrentSides(ll);
                else
                    LegLocations(CurrentLegs(ll),1) = kk;
                    LegSides(CurrentLegs(ll),1) = CurrentSides(ll);
                end
            else
                LegLocationsNeg(-CurrentLegs(ll),2) = kk;
                LegSidesNeg(-CurrentLegs(ll),2) = CurrentSides(ll);
            end
        end
    end
    
    if any(any(LegLocations==0))
        error('Affirmation Error: We should have found 2 ends for all Legs in Loop Part');
    end
    
    LegEndLocations = cell(size(LegLocations));
    
    if isempty(UniqueLegsNeg)
        FinishedTensors = false([TotalTensors(end),1]);
        LockedTensors = false([TotalTensors(end),1]);
        TensorPath = cell([TotalTensors(end),1]);
    else
        FinishedTensors = false([TotalTensors(end)+1,1]);
        LockedTensors = false([TotalTensors(end)+1,1]);
        TensorPath = cell([TotalTensors(end)+1,1]);
    end
    LoopCounter = 0;
    
    PathOff = zeros([0,2]);
    PathOffCounter = [];
    Counter = 0;
    UsingTensors = zeros([1,0]);%unused tensors
    while ~all(FinishedTensors) %if we haven't finished all tensors.
        
        %details regarding Pathoff, it is a list of leg locations in a list
        %of tensors, the first being the tensor number and the second being
        %the leg number in this ordering (after we use RotateTensors)
        %if we are refering to tensor locations then the LegNumber is set
        %to zero.
        %multiple rows refers to lists coming off other legs
        
        
        SeedLegLocation = cell([1,0]);
        SeedTensors = zeros([1,0]);
        
        %search through last entry for possible layers that we can create
        %and connect any allowed connections
        
        if ~isempty(UsingTensors)
            %first check that there are still free legs
            Numbers = 1:size(CyclicLegsBraidedLabelsUser{UsingTensors(end)},2);
            FreeLegs = CyclicLegsBraidedLabelsUser{UsingTensors(end)}(~CyclicLegsBraidedUsed{UsingTensors(end)});
            Numbers = Numbers(~CyclicLegsBraidedUsed{UsingTensors(end)});
            CurrentPath = TensorPath{UsingTensors(end)};
            CurrentTensor = UsingTensors(end);
            
            
            if isempty(FreeLegs)
                FinishedTensors(UsingTensors(end)) = true;
                UsingTensors(end) = [];
                
            else
                
                CountFreeLegs = sum(LegLocations(FreeLegs,:)'==CurrentTensor,1);
                
                %first count all the traces
                TraceLegs = unique(FreeLegs(CountFreeLegs == 2));
                FreeLegs(CountFreeLegs == 2) = [];
                TraceLegsLoc = Numbers(CountFreeLegs == 2);
                Numbers(CountFreeLegs == 2) = [];
                
                %Now set the LegEndLocations for the traced legs:
                
                CyclicLegsBraidedUsed{UsingTensors(end)}(TraceLegsLoc) = true;
                CyclicLegsBraidedFirstSide{UsingTensors(end)}(TraceLegsLoc) = true;
                if ~isempty(TraceLegs)
                    for tt = TraceLegs
                        LegLoc = find(CyclicLegsBraidedLabelsUser{UsingTensors(end)} == tt);
                        if numel(LegLoc)~=2
                             error('Affirmation Error: We didn''t find the right number of Legs');
                        end
                        
                        LegEndLocations{tt,1} = CurrentPath;
                        LegEndLocations{tt,1}(end,2) = LegLoc(1);
                        LegEndLocations{tt,2} = CurrentPath;
                        LegEndLocations{tt,2}(end,2) = LegLoc(2);
                    end
                end
                
                OtherTensors = sum(LegLocations(FreeLegs,:).*(LegLocations(FreeLegs,:)~=UsingTensors(end)),2)';
                
                LegsLockedToSomewhere = LockedTensors(OtherTensors);
                CyclicLegsBraidedUsed{UsingTensors(end)}(Numbers) = true;
                %CyclicLegsBraidedFirstSide{UsingTensors(end)}(Numbers(~LegsLockedToSomewhere)) = true;
                
                SeedTensors = OtherTensors(~LegsLockedToSomewhere);
                SeedLegLocation = repmat(TensorPath(UsingTensors(end)),[1,sum(~LegsLockedToSomewhere)]);
                SeedLegs = Numbers(~LegsLockedToSomewhere);
                SeedLegsLabels = FreeLegs(~LegsLockedToSomewhere);
                
                for kk = 1:numel(SeedLegLocation);
                    SeedLegLocation{kk}(end,2) = SeedLegs(kk);
                    %spawn from a certain Leg
                end
                
                FreeLegs(~LegsLockedToSomewhere) = [];
                Numbers(~LegsLockedToSomewhere) = [];
                OtherTensors(~LegsLockedToSomewhere) = [];
                
                %now we have the details for the next layers we are going
                %to generate, specifically built out of:
                % - SeedLegLocations
                % - SeedTensors
                
                %now to put all the details in the fixed Locations, as we
                %are going from left to right this means if they haven't
                %been set yet the current tensor is left-most (and the
                %first entry of LegEndLocations
                
                
                if sum(LegsLockedToSomewhere)>0
                    for aa = 1:numel(Numbers)
                        tt = FreeLegs(aa);
                        LegEndLocations{tt,1} = CurrentPath;
                        LegEndLocations{tt,1}(end,2) = Numbers(aa); %we have assumed that this is the smallest one
                        
                        %now other point
                        LegEndLocations{tt,2} = TensorPath{OtherTensors(aa)};
                        Temp = find(CyclicLegsBraidedLabelsUser{OtherTensors(aa)} == tt);
                        if numel(Temp)~=1
                            error('Affirmation Error: Wrong number of times found for other term');
                        end
                        LegEndLocations{tt,2}(end,2) = Temp;
                        
                        if GTPath(LegEndLocations{tt,2}, LegEndLocations{tt,1})
                            CyclicLegsBraidedFirstSide{CurrentTensor}(Numbers(aa)) = true;
                        else
                            CyclicLegsBraidedFirstSide{OtherTensors(aa)}(Temp) = true;
                            LegEndLocations(tt,1:2) = LegEndLocations(tt,[2,1]);
                        end
                    end
                end
                
                %now remove UsedTensors(end)
                CyclicLegsBraidedUsed{CurrentTensor} = true(size(CyclicLegsBraidedUsed{CurrentTensor}));
                FinishedTensors(CurrentTensor) = true;
                UsingTensors(end) = [];
            end
        else     
            %if it is empty then set a new tensor to be the source of the list
            if ~isequal(LockedTensors, FinishedTensors)
                error('Affirmation Error: we should have all tensors which are locked are finished when we have no Using tensors left.')
            end
            SeedTensors = find(~LockedTensors,1,'first');
            if numel(SeedTensors) ~=1
                error('Affirmation Error: There should be tensors to find which aren''t locked')
            end
            Counter = Counter+1;
            CurrentPath = [Counter,0];
            SeedLegLocation = {CurrentPath}; %note that this zero means that we aren't coming off a leg.
            CurrentTensor = 0;
        end
        
        
        
        %take the list of new layer seeds, note that if the seed gets used
        %previously (during a previous layer generated from a seed then we
        %skip it
        
        if ~isempty(SeedTensors)
            %first check if we have that the current Tensor is zero (in
            %which case we are just adding a new bottom layer)
            
            AddToCurrentUsing = zeros([1,0]);
            
            if CurrentTensor == 0
                
                if numel(SeedTensors) ~= 1
                    error('Affirmation Error: We have the wrong number of Seed Tensors for a new bottom layer')
                end
                
                RotateTensor(SeedTensors) = 1;
                
                if ~all(CyclicLegsBraidedUsed{SeedTensors})
                
                    LegLoc = find(~CyclicLegsBraidedUsed{SeedTensors},1,'last');
                    if numel(LegLoc) ~= 1
                        error('Affirmation Error: can''t find a last leg to create the next disjoint set')
                    end
                    
                    %CyclicLegsBraidedUsed{SeedTensors}(LegLoc) = true;
                    %CyclicLegsBraidedFirstSide{SeedTensors}(LegLoc) = true;
                    
                    %The following are kept the same
                    % - CyclicLegsBraidedExternalLeg(SeedTensors)
                    % - CyclicLegsBraidedLabelsUser(SeedTensors)
                    
                    %Now work out which other tensor is connected to this 
                    
                    CountNewLayer = 1;
                    
                    FlagContinue = true;
                    WorkingTensors = SeedTensors;
                    LockedTensors(SeedTensors) = true;
                    WorkingExternalLeg = CyclicLegsBraidedExternalLeg{WorkingTensors};
                    WorkingLegs = CyclicLegsBraidedLabelsUser{WorkingTensors};
                    AddToCurrentUsing = [AddToCurrentUsing,WorkingTensors];
                    while FlagContinue
                        TensorOptions = LegLocations(CyclicLegsBraidedLabelsUser{WorkingTensors}(LegLoc),:);
                        if ~all(TensorOptions == WorkingTensors)
                            if ~any(TensorOptions == WorkingTensors)
                                error('Affirmation Errors: Something went wrong with LegLocations and SeedTensors')
                            end
                            NextTensor = TensorOptions(TensorOptions~=WorkingTensors);
                            NextLeg = CyclicLegsBraidedLabelsUser{WorkingTensors}(LegLoc);
                            
                            CheckTensors = WorkingTensors;
                            CheckLeg = LegLoc;
                            
                            CyclicLegsBraidedUsed{WorkingTensors}(LegLoc) = true;
                            TensorPath{WorkingTensors} = [CurrentPath;[CountNewLayer,0]];
                            LegEndLocations{NextLeg, 1} = [CurrentPath;[CountNewLayer,LegLoc]];
                            CountNewLayer = CountNewLayer+1;
                            FlagContinue = false;
                            
                        else %then this is a trace
                            LegLocTemp = find(WorkingLegs == WorkingLegs(LegLoc));
                            if isempty(LegLocTemp)
                                error('Affirmation Error: This doesn''t make sense, the LegLocTemp shouldn''t be empty');
                            end
                            
                            LegLoc = find(~WorkingExternalLeg(1:(LegLocTemp-1)),1,'last');
                            
                            if numel(LegLocTemp)==0
                                NextTensor = [];
                                FlagContinue = false;
                            end
                                
                        end
                    end
                    
                    while ~isempty(NextTensor)
                        if LockedTensors(NextTensor)
                            LegEndLocations{NextLeg, 2} = TensorPath{NextTensor};
                            
                            LegLoc = find(CyclicLegsBraidedLabelsUser{NextTensor} == NextLeg);
                            if numel(LegLoc) ~=1
                                error('Affirmation Error: This is the wrong number of entries for rotation')
                            end
                            LegEndLocations{NextLeg, 2}(end,2) = LegLoc;
                            
                            %2 is the locked tensor, so to check we want 2
                            %to be greater then (After) 1, if that is
                            %alredy true then we need to set the checked
                            %locations to true.
                            
                            if GTPath(LegEndLocations{NextLeg, 2},LegEndLocations{NextLeg, 1})
                                CyclicLegsBraidedFirstSide{CheckTensors}(CheckLeg) = true;
                            else
                                CyclicLegsBraidedFirstSide{NextTensor}(LegLoc) = true;
                                LegEndLocations(NextLeg, 1:2) = LegEndLocations(NextLeg, [2,1]);
                            end
                            
                            CyclicLegsBraidedUsed{NextTensor}(LegLoc) = true;
                            break;
                        end
                        
                        FlagContinue = true;
                        WorkingTensors = NextTensor;
                        
                        RotateTensorTemp = find(CyclicLegsBraidedLabelsUser{WorkingTensors} == NextLeg);
                        if numel(RotateTensorTemp) ~=1
                            error('Affirmation Error: This is the wrong number of entries for rotation')
                        end
                        
                        CyclicLegsBraidedLabelsUser{WorkingTensors} = CyclicLegsBraidedLabelsUser{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        CyclicLegsBraidedFirstSide{WorkingTensors} = CyclicLegsBraidedFirstSide{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        CyclicLegsBraidedExternalLeg{WorkingTensors} = CyclicLegsBraidedExternalLeg{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        CyclicLegsBraidedUsed{WorkingTensors} = CyclicLegsBraidedUsed{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        RotateTensor(WorkingTensors) = RotateTensorTemp;
                        TensorPath{WorkingTensors} = [CurrentPath;[CountNewLayer,0]];
                        LegEndLocations{NextLeg, 2} = [CurrentPath;[CountNewLayer,1]];%if this is new it is on the right
                        CyclicLegsBraidedFirstSide{CheckTensors}(CheckLeg) = true; %correct the first side
                        
                        CyclicLegsBraidedUsed{WorkingTensors}(1) = true;
                        WorkingExternalLeg = CyclicLegsBraidedExternalLeg{WorkingTensors};
                        WorkingLegs = CyclicLegsBraidedLabelsUser{WorkingTensors};
                        LegLoc = find(~WorkingExternalLeg,1,'last');
                        AddToCurrentUsing = [AddToCurrentUsing,WorkingTensors];
                        LockedTensors(WorkingTensors) = true;
                        if isempty(LegLoc)
                            error('Affimation Error: We should at least find the initial leg')
                        end
                        if LegLoc == 1
                            FlagContinue = false;
                            NextTensor = [];
                            %there are no more tensors to add to this.
                        end
                        while FlagContinue
                            TensorOptions = LegLocations(CyclicLegsBraidedLabelsUser{WorkingTensors}(LegLoc),:);
                            if ~all(TensorOptions == WorkingTensors)
                                if ~any(TensorOptions == WorkingTensors)
                                    error('Affirmation Errors: Something went wrong with LegLocations and SeedTensors')
                                end
                                NextTensor = TensorOptions(TensorOptions~=WorkingTensors);
                                NextLeg = CyclicLegsBraidedLabelsUser{WorkingTensors}(LegLoc);
                                CyclicLegsBraidedUsed{WorkingTensors}(LegLoc) = true;
                                TensorPath{WorkingTensors} = [CurrentPath;[CountNewLayer,0]];
                                LegEndLocations{NextLeg, 1} = [CurrentPath;[CountNewLayer,LegLoc]];
                                CheckTensors = WorkingTensors;
                                CheckLeg = LegLoc;
                                CountNewLayer = CountNewLayer+1;
                                FlagContinue = false;
                                
                            else %then this is a trace
                                LegLocTemp = find(WorkingLegs == WorkingLegs(LegLoc));
                                if isEmpty(LegLocTemp)
                                    error('Affirmation Error: This doesn''t make sense, the LegLocTemp shouldn''t be empty');
                                end
                                
                                LegLoc = find(~WorkingExternalLeg(1:(LegLocTemp-1)),1,'last');
                                
                                if numel(LegLocTemp)==0
                                    NextTensor = [];
                                    FlagContinue = false;
                                end
                                if LegLocTemp==1
                                    NextTensor = [];
                                    FlagContinue = false;
                                end
                                    
                            end
                        end
                    end
                else
                    if any(CyclicLegsBraidedLabelsUser{WorkingTensors}>0)
                        error('Affirmation Error: This tensor should be completely disconnected from all other tensors (there should be no internal legs)')
                    end
                    TensorPath{SeedTensors} = [CurrentPath;[1,0]];
                    FinishedTensors(SeedTensors) = true;
                    LockedTensors(SeedTensors) = true;
                end
                
            else
                CheckTensorsRaw = CurrentTensor;
                
                for nn = 1:length(SeedTensors)
                    CheckTensors = CheckTensorsRaw;
                    CheckLeg = SeedLegs(nn);
                    if LockedTensors(SeedTensors(nn))
                        continue;
                    end
                    %CurrentPath = SeedTensors{NewLayerSeeds(nn)};
                    % - SeedLegLocation
                    % - SeedLegLabels
                    % - SeedTensors
                    
                    CountNewLayer = 1;
                    NextTensor = SeedTensors(nn);
                    NextLeg = SeedLegsLabels(nn);
                    CurrentPath = SeedLegLocation{nn};
                    LegEndLocations{NextLeg,1} = CurrentPath;
                    
                    
                    while ~isempty(NextTensor)
                        if LockedTensors(NextTensor)
                            LegEndLocations{NextLeg, 2} = TensorPath{NextTensor};
                            
                            LegLoc = find(CyclicLegsBraidedLabelsUser{NextTensor} == NextLeg);
                            if numel(LegLoc) ~=1
                                error('Affirmation Error: This is the wrong number of entries for rotation')
                            end
                            LegEndLocations{NextLeg, 2}(end,2) = LegLoc;
                            
                            %2 is the locked tensor, so to check we want 2
                            %to be greater then (After) 1, if that is
                            %alredy true then we need to set the checked
                            %locations to true.
                            
                            if GTPath(LegEndLocations{NextLeg, 2},LegEndLocations{NextLeg, 1})
                                CyclicLegsBraidedFirstSide{CheckTensors}(CheckLeg) = true;
                            else
                                CyclicLegsBraidedFirstSide{NextTensor}(LegLoc) = true;
                                LegEndLocations(NextLeg, 1:2) = LegEndLocations(NextLeg, [2,1]);
                            end
                            
                            CyclicLegsBraidedUsed{NextTensor}(LegLoc) = true;
                            break;
                        end
                        
                        FlagContinue = true;
                        WorkingTensors = NextTensor;
                        
                        RotateTensorTemp = find(CyclicLegsBraidedLabelsUser{WorkingTensors} == NextLeg);
                        if numel(RotateTensorTemp) ~=1
                            error('Affirmation Error: This is the wrong number of entries for rotation')
                        end
                        
                        CyclicLegsBraidedLabelsUser{WorkingTensors} = CyclicLegsBraidedLabelsUser{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        CyclicLegsBraidedFirstSide{WorkingTensors} = CyclicLegsBraidedFirstSide{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        CyclicLegsBraidedExternalLeg{WorkingTensors} = CyclicLegsBraidedExternalLeg{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        CyclicLegsBraidedUsed{WorkingTensors} = CyclicLegsBraidedUsed{WorkingTensors}([RotateTensorTemp:end,1:(RotateTensorTemp-1)]);
                        RotateTensor(WorkingTensors) = RotateTensorTemp;
                        TensorPath{WorkingTensors} = [CurrentPath;[CountNewLayer,0]];
                        LegEndLocations{NextLeg, 2} = [CurrentPath;[CountNewLayer,1]];
                        CyclicLegsBraidedFirstSide{CheckTensors}(CheckLeg) = true; %correct the first side
                        
                        CyclicLegsBraidedUsed{WorkingTensors}(1) = true;
                        WorkingExternalLeg = CyclicLegsBraidedExternalLeg{WorkingTensors};
                        WorkingLegs = CyclicLegsBraidedLabelsUser{WorkingTensors};
                        LegLoc = find(~WorkingExternalLeg,1,'last');
                        AddToCurrentUsing = [AddToCurrentUsing,WorkingTensors];
                        LockedTensors(WorkingTensors) = true;
                        if isempty(LegLoc)
                            error('Affimation Error: We should at least find the initial leg')
                        end
                        if LegLoc == 1
                            FlagContinue = false;
                            NextTensor = [];
                            %there are no more tensors to add to this.
                        end
                        while FlagContinue
                            TensorOptions = LegLocations(CyclicLegsBraidedLabelsUser{WorkingTensors}(LegLoc),:);
                            if ~all(TensorOptions == WorkingTensors)
                                if ~any(TensorOptions == WorkingTensors)
                                    error('Affirmation Errors: Something went wrong with LegLocations and SeedTensors')
                                end
                                NextTensor = TensorOptions(TensorOptions~=WorkingTensors);
                                NextLeg = CyclicLegsBraidedLabelsUser{WorkingTensors}(LegLoc);
                                
                                CheckTensors = WorkingTensors;
                                CheckLeg = LegLoc;
                                
                                CyclicLegsBraidedUsed{WorkingTensors}(LegLoc) = true;
                                TensorPath{WorkingTensors} = [CurrentPath;[CountNewLayer,0]];
                                LegEndLocations{NextLeg, 1} = [CurrentPath;[CountNewLayer,LegLoc]];
                                CountNewLayer = CountNewLayer+1;
                                FlagContinue = false;
                                
                            else %then this is a trace
                                LegLocTemp = find(WorkingLegs == WorkingLegs(LegLoc));
                                if isEmpty(LegLocTemp)
                                    error('Affirmation Error: This doesn''t make sense, the LegLocTemp shouldn''t be empty');
                                end
                                
                                LegLoc = find(~WorkingExternalLeg(1:(LegLocTemp-1)),1,'last');
                                
                                if numel(LegLocTemp)==0
                                    NextTensor = [];
                                    FlagContinue = false;
                                end
                                if LegLocTemp==1
                                    NextTensor = [];
                                    FlagContinue = false;
                                end
                                    
                            end
                        end
                    end
                    
                end
                
            end
            
            UsingTensors = [UsingTensors,AddToCurrentUsing(end:-1:1)];
            
        end
        
    end
    
    %here we are going to work out the Type of each leg
    
    LegType = cell(size(CyclicLegsBraidedTopBottomSide));
    
    for kk = 1:numel(CyclicLegsBraidedTopBottomSide)
        Numbers = 1:numel(CyclicLegsBraidedTopBottomSide{kk});
        TopNumbers = Numbers(CyclicLegsBraidedTopBottomSide{kk} == +1);
        BottomNumbers = Numbers(CyclicLegsBraidedTopBottomSide{kk} == -1);
        %BottomNumbers = BottomNumbers(end:-1:1);
        LegType{kk} = zeros([1,numel(Numbers)]);
                
        if RotateTensor(kk) == 1
            LegType{kk}(TopNumbers) = 1; % Good
            LegType{kk}(BottomNumbers) = 4; %UpRight
            
            %LegsOutOrder = [TopNumbers, BottomNumbers];
        elseif RotateTensor(kk)>numel(TopNumbers)
            %then we don't do any rotations on the top
            LegType{kk}(TopNumbers) = 1; %Good
            LegType{kk}(BottomNumbers(1:(RotateTensor(kk)-numel(TopNumbers)-1))) = 4; %UpRight
            
            LegType{kk}(BottomNumbers((RotateTensor(kk)-numel(TopNumbers)):end)) = 3; %UpLeft
            
            %LegsOutOrder = [BottomNumbers((RotateTensor(kk)-numel(TopNumbers)):end), TopNumbers, BottomNumbers(1:(RotateTensor(kk)-numel(TopNumbers)-1))];
        else
            %then we have all bottoms are upLeft
            LegType{kk}(BottomNumbers) = 3; %UpLeft
            LegType{kk}(TopNumbers(RotateTensor(kk):end)) = 2; %Rotate
            
            LegType{kk}(TopNumbers(1:(RotateTensor(kk)-1))) = 1; %Good
            
            %LegsOutOrder = [TopNumbers(RotateTensor(kk):end), BottomNumbers, TopNumbers(1:(RotateTensor(kk)-1))];
        end
        
    end
    
    if isempty(UniqueLegsNeg)
        ModificationsSinglePhase = cell([numel(CyclicLegsBraidedLabelsUserOld),1]);
        ModificationsEnvironmentsPhase = cell([numel(CyclicLegsBraidedLabelsUserOld),1]);
    else
        ModificationsSinglePhase = cell([numel(CyclicLegsBraidedLabelsUserOld)-1,1]);
        ModificationsEnvironmentsPhase = cell([numel(CyclicLegsBraidedLabelsUserOld)-1,1]);
    end
    
    for kk = 1:numel(ModificationsSinglePhase)
        ModificationsSinglePhase{kk} = cell([0,2]);
        ModificationsEnvironmentsPhase{kk} = zeros([0,2]);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now we know the planar diagram
    %now check that the loops give rise to planar diagrams
    %do this by comparing all legs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for kk = 1:size(LegEndLocations,1)
        if GTPath(LegEndLocations{kk,1}, LegEndLocations{kk,2})
            %if 1 is greater (more right) then 2
            error('Affirmation Error: I should have the legs ordered so that the first end is the more left most end');
        end
    end
    
    for kk = 1:size(LegEndLocations,1)
        
        for ll = (kk+1):size(LegEndLocations,1)
            
            Direction11 = GTPath(LegEndLocations{kk,1}, LegEndLocations{ll,1});
            Direction12 = GTPath(LegEndLocations{kk,1}, LegEndLocations{ll,2});
            Direction22 = GTPath(LegEndLocations{kk,2}, LegEndLocations{ll,2});
            Direction21 = GTPath(LegEndLocations{kk,2}, LegEndLocations{ll,1});
            
            if ~((~Direction12) || Direction21)
                if Direction11==Direction22
                    error('Error: This is not a planar diagram, I cannot contract this correctly')
                end
            end
        end
        
    end
    
    %now make sure LegSides and LegLocations are ordered the same as
    %LegEndLocations
    
    for kk = 1:size(LegEndLocations,1)
        LegEndTemp = LegEndLocations{kk,1};
        LegEndTemp(end,2) = 0;
        if ~isequal(LegEndTemp,TensorPath{LegLocations(kk,1)})
            LegSides(kk,:) = LegSides(kk,[2,1]);
            LegLocations(kk,:) = LegLocations(kk,[2,1]);
            if ~isequal(LegEndTemp,TensorPath{LegLocations(kk,1)})
                error('Affirmation Error: The LegLocations are weird')
            end
        end
    end
    
    LegTypeLocations = zeros(size(LegLocations));
    for kk = 1:size(LegTypeLocations,1)
        if LegLocations(kk,1)~=LegLocations(kk,2)
            LocA = LegLocations(kk,1);
            TypeA = LegType{LocA}(CyclicLegsBraidedLabelsUserOld{LocA}==kk);
            LocB = LegLocations(kk,2);
            TypeB = LegType{LocB}(CyclicLegsBraidedLabelsUserOld{LocB}==kk);
            
            %A is the first tensor, then B is the second tensor
            
            if numel(TypeA)~=1
                error('Affirmation Error: TypeA is not a single value')
            end
            if numel(TypeB)~=1
                error('Affirmation Error: TypeB is not a single value')
            end
            
            if sum([TypeA==[3,4],TypeB==[3,4]])~=1
                error('Affirmation Error: They are both on the same side')
            end
            
            LegTypeLocations(kk,:) = [TypeA,TypeB];
            %then need to modify
            if ~isempty(UniqueLegsNeg)&&LocB==numel(CyclicLegsBraidedLabelsUserOld)
                %special case: B is the output tensor, so therefor
                %apply to A
                
                N = find(CyclicLegsBraidedLabelsUser{LocA}==kk);
                Dir = ChargeDirectionsTensorDetails{LocA}(StructureTensorDetails{LocA}==-CyclicLegsBraidedTensorDetails{LocA}(N));
                
                if TypeA==4
                    ModificationsSinglePhase{LocA} = [ModificationsSinglePhase{LocA};{'pDown',[-N;Dir]}];
                    ModificationsEnvironmentsPhase{LocA} = [ModificationsEnvironmentsPhase{LocA};kk,Dir];
                elseif TypeB == 3
                    ModificationsSinglePhase{LocA} = [ModificationsSinglePhase{LocA};{'pDown',[-N;Dir]}];
                    ModificationsEnvironmentsPhase{LocA} = [ModificationsEnvironmentsPhase{LocA};kk,Dir];
                end
                %if TypeB == 4 or TypeA == 3 then get the identity
            else
                
                N = find(CyclicLegsBraidedLabelsUser{LocB}==kk);
                Dir = ChargeDirectionsTensorDetails{LocB}(StructureTensorDetails{LocB}==-CyclicLegsBraidedTensorDetails{LocB}(N));
                
                if TypeA==4
                    ModificationsSinglePhase{LocB} = [ModificationsSinglePhase{LocB};{'pDown',[-N;-Dir]}];
                    ModificationsEnvironmentsPhase{LocB} = [ModificationsEnvironmentsPhase{LocB};kk,-Dir];
                elseif TypeB == 3
                    ModificationsSinglePhase{LocB} = [ModificationsSinglePhase{LocB};{'pDown',[-N;-Dir]}];
                    ModificationsEnvironmentsPhase{LocB} = [ModificationsEnvironmentsPhase{LocB};kk,-Dir];
                end
                %if TypeB == 4 or TypeA == 3 then get the identity
            end
            
        end
    end
    
    %now add modififers based on this:
    
    
    
    
    
    
    %add a check in here to make sure that first agrees with
    %LegEndLocations
    
    for kk = 1:numel(CyclicLegsBraidedFirstSide)
        for aa = 1:numel(CyclicLegsBraidedFirstSide{kk})
            if CyclicLegsBraidedLabelsUser{kk}(aa)>0
               ThisTensor = LegLocations(CyclicLegsBraidedLabelsUser{kk}(aa),:)==kk;
               %check that if the location is first then FirstSide must be
               %true
               if ThisTensor(1)&&~ThisTensor(2)&&~CyclicLegsBraidedFirstSide{kk}(aa)
                   error('Affirmation Error: There is a problem with CyclicLegsBraidedFirstSide, it doesn''t agree with LegLocations')
               end
               %check that if the location is second then FirstSide must be
               %false
               if ThisTensor(1)&&~ThisTensor(2)&&~CyclicLegsBraidedFirstSide{kk}(aa)
                   error('Affirmation Error: There is a problem with CyclicLegsBraidedFirstSide, it doesn''t agree with LegLocations')
               end
                
            end
        end
    end
    
    
    if ~isempty(UniqueLegsNeg)
        RotateTensorInitial = RotateTensor(1:(end-1));
        NegativeRotateTensor = RotateTensor(end);
        CyclicLegsBraidedFirstSide = CyclicLegsBraidedFirstSide(1:(end-1));
        
        for kk = 1:NumberTensors
            CyclicLegsBraidedExternalLeg{kk} = CyclicLegsBraidedExternalLegNew{kk}([RotateTensor(kk):end,1:(RotateTensor(kk)-1)]);
        end
    else
        RotateTensorInitial = RotateTensor;
        NegativeRotateTensor = 0;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now we have the structure the next thing is to rotate the tensors so
    %that this is true (this is the restructuring)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %first work out if we need rotation modifiers
    
    RotationModifiersInitial = SymTensor.FromStandard(Tensors, RotateTensorInitial, false);
    
    for kk = 1:numel(RotationModifiersInitial)
        RotationModifiersInitial{kk} = [RotationModifiersInitial{kk};ModificationsSinglePhase{kk}];
    end
    
    %now consider if we need to add phases at the begining to account for
    %putting in the caps
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now do the Contraction corrections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %in this part I assume I can use:
    % X- ActionTensorsLocations
    % - ActionTypeTotal
    % X- ActionTensorsPermute
    % X- ActionTensorsConvert
    % - ActionNumberContractedLegs
    % X- C_Structure
    % X- C_ChargeDirections
    % X- DropTensor
    % X- OldTrueLocations
    % X- OldTrueConvert
    % X- OldTruePermute
    % - OldTrueChargeSides
    % - OldTrueStructure
    % - OldTrueChargeDirections
    % X- OldTrueStoredLocations
    %as well as:
    % - ActionTotalTensors
    % - TotalTensors
    % - ActualLegsTotal
    % - OldTrue
    % - LegOrder
    % - ActionTensorsStructure
    % - ActionTensorsChargeDirections

    
    
    %in this part I am updating:
    % - ActionTensorsLocations
    % - ActionTensorsPermute
    % - ActionTensorsConvert
    % - C_Structure
    % - C_ChargeDirections
    % - OldTrueConvert
    % - OldTruePermute
    % - OldTrueChargeSides
    % - OldTrueStructure
    % - OldTrueChargeDirections
    % and also:
    % - ActionTotalTensors
    % - ActualLegsTotal*
    % - ActionTensorsStructure*
    % - ActionTensorsChargeDirections*
    
    %first reset inital parameters
    
    ActionTensorsStructure = cell(size(ActionTotalTensors));
    ActionTensorsChargeDirections = cell(size(ActionTotalTensors));
    ActionTensorsChargeDirectionsExternal = cell(size(ActionTotalTensors));
    ActionTensorsChargeDirectionsExternalOriginal = cell(size(ActionTotalTensors));
    ActualLegsTotal = cell(size(ActionTotalTensors));
    TensorPathTotal = cell(size(ActionTotalTensors));
    
    
    ActionCyclicLegsConnectedTensor = cell(size(ActionTotalTensors));
    ActionCyclicLegsIsTrace = cell(size(ActionTotalTensors));
    ActionCyclicLegsIsTraceShaded = cell(size(ActionTotalTensors));
    ActionCyclicLegsBraidedFirstSide = cell(size(ActionTotalTensors));
    ActionCyclicLegsBraidedCountFlipSide = cell(size(ActionTotalTensors));
    ActionCyclicLegsBraidedExternalLeg = cell(size(ActionTotalTensors));
    ActionTensorsChargeDirectionsSides = cell(size(ActionTotalTensors));
    ActionTensorsChargeDirections = cell(size(ActionTotalTensors));
    
    ActionTensorsLegLocations = cell(size(ActionTotalTensors));
    ActionTensorsLegPosModified = cell(size(ActionTotalTensors));
    ActionTensorsClumptedTensors = cell(size(ActionTotalTensors));
    ActionTensorsLegEndLocations = cell(size(ActionTotalTensors));
    ActionTensorsLegEndLocationsOriginal = LegEndLocations;
    ActionTensorsDictionary = repmat({cell([0,2])},size(ActionTotalTensors));
    
    CyclicLegsBraidedLabelsUserFull = CyclicLegsBraidedLabelsUser;
    
    NewOutTensorSpherical = 0;
    if ~isempty(RenameNeg)
        OutTensorSpherical = numel(CyclicLegsBraidedLabelsUser);
        
        LegLocations(LegLocations == OutTensorSpherical) = NewOutTensorSpherical;
        TensorSphericalPath = TensorPath{OutTensorSpherical};
        TensorSphericalNegsLegs = CyclicLegsBraidedLabelsUser{OutTensorSpherical};
        for kk = 1:numel(TensorSphericalNegsLegs)
            TensorSphericalNegsLegs(kk) = RenameNeg(1,RenameNeg(2,:)==TensorSphericalNegsLegs(kk));
        end
        %TensorPath(max(max(ActionTotalTensors))+1) = TensorPath(OutTensorSpherical);
        
    end
    for nn = 1:NumberTensors
        if ~isempty(UniqueLegsNeg)
            for kk = 1:size(RenameNeg,2)
                CyclicLegsBraidedLabelsUser{nn}(CyclicLegsBraidedLabelsUser{nn}==RenameNeg(2,kk)) = RenameNeg(1,kk);
            end
        end
        
        Temp = ActionTotalTensors==nn;
        
        
        ActionTensorsLegLocations(Temp) = {LegLocations};
        ActionTensorsLegPosModified(Temp) = {LegLocations == nn};
        ActionTensorsLegEndLocations(Temp) = {LegEndLocations};
        ActionTensorsClumptedTensors(Temp) = {nn};
        
        ActualLegsTotal(Temp) = CyclicLegsBraidedLabelsUser(nn);
        TensorPathTotal(Temp) = TensorPath(nn);
        
        
        
        
        
        NumberLegs = numel(CyclicLegsBraidedLabelsUser{nn});
        [~,Index] = sort(StructureTensorDetails{nn}(StructureTensorDetails{nn}<0),'descend');
        ChargeDirectionsExternal = ChargeDirectionsTensorDetails{nn}(StructureTensorDetails{nn}<0);
        
        ChargeDirectionsExternal = ChargeDirectionsExternal(Index);
        ChargeDirectionsExternalFlat = ChargeDirectionsExternal(CyclicLegsBraidedTensorDetails{nn}([RotateTensor(nn):end,1:(RotateTensor(nn)-1)]));
        
        ChargeSides = ChargeSideTensorDetails{nn}(CyclicLegsBraidedTensorDetails{nn}([RotateTensor(nn):end,1:(RotateTensor(nn)-1)]));
        
        %IndexFrom
        %IndexFromOriginal = Index(CyclicLegsBraidedTensorDetails{nn}([RotateTensor(nn):end,1:(RotateTensor(nn)-1)]));
        
        %ChargeDirectionsExternal = ChargeDirectionsExternal(IndexFromOriginal);
        %ChargeSides = ChargeSideTensorDetails{nn}(IndexFromOriginal);
        
        
        
        
        
        if NumberLegs < 2
            StructureTemp = [];
            ChargeDirectionsTemp = [];
        elseif NumberLegs == 2
            StructureTemp = [-1;-2];
            ChargeDirectionsTemp = ChargeDirectionsExternal(CyclicLegsBraidedTensorDetails{nn}(mod([0;1]+RotateTensor(nn),NumberLegs)+1));
        else
            StructureTemp = [[NumberLegs-2;-NumberLegs],[-1;-2],[1:(NumberLegs-3);-(3:(NumberLegs-1))]];
            ChargeDirectionsTemp = [[+1;ChargeDirectionsExternal(CyclicLegsBraidedTensorDetails{nn}(mod(NumberLegs-3+RotateTensor(nn),NumberLegs)+1))],...
                ChargeDirectionsExternal(CyclicLegsBraidedTensorDetails{nn}(mod([0;1]+RotateTensor(nn),NumberLegs)+1)),...
                [ones([1,NumberLegs-3]);reshape(ChargeDirectionsExternal(CyclicLegsBraidedTensorDetails{nn}(mod((2:(NumberLegs-2))+RotateTensor(nn),NumberLegs)+1)),[1,NumberLegs-3])]];
        end
        
        ActionTensorsChargeDirectionsExternal(Temp) = repmat({reshape(ChargeDirectionsExternalFlat,[1,numel(ChargeDirectionsExternal)])}, [sum(sum(Temp)),1]);
        ActionTensorsChargeDirectionsExternalOriginal(Temp) = repmat({reshape(ChargeDirectionsExternal,[1,numel(ChargeDirectionsExternal)])}, [sum(sum(Temp)),1]);
        ActionTensorsStructure(Temp) = repmat({StructureTemp}, [sum(sum(Temp)),1]);
        ActionTensorsChargeDirections(Temp) = repmat({ChargeDirectionsTemp}, [sum(sum(Temp)),1]);
        
        ActionCyclicLegsBraidedFirstSide(Temp) = CyclicLegsBraidedFirstSide(nn);
        ActionCyclicLegsBraidedCountFlipSide(Temp) = {zeros(size(CyclicLegsBraidedFirstSide{nn}))};
        ActionCyclicLegsBraidedExternalLeg(Temp) = CyclicLegsBraidedExternalLeg(nn);
        
        OtherTensor = LegLocations(CyclicLegsBraidedLabelsUser{nn}(~CyclicLegsBraidedExternalLeg{nn}),:)';
        TraceLegs = OtherTensor(1,:)==OtherTensor(2,:);
        if any(any(OtherTensor(:,TraceLegs)~=nn))
            error('Affimation Error: something wrong with working out traces')
        end
        OtherTensor(1,~TraceLegs) = OtherTensor(OtherTensor~=nn);
        OtherTensor = OtherTensor(1,:);
        
        IsTrace = false(size(CyclicLegsBraidedExternalLeg{nn}));
        ConnectedTensor = zeros(size(CyclicLegsBraidedExternalLeg{nn}));
        
        IsTrace(~CyclicLegsBraidedExternalLeg{nn}) = TraceLegs;
        ConnectedTensor(~CyclicLegsBraidedExternalLeg{nn}) = OtherTensor;
        
        TraceNumbers = 1:length(IsTrace);
        TraceNumbers = TraceNumbers(IsTrace);
        IsTraceShaded = false(size(IsTrace));
        for kk = TraceNumbers
            EndLocs = find(CyclicLegsBraidedLabelsUser{nn} == CyclicLegsBraidedLabelsUser{nn}(kk));
            if numel(EndLocs) ~=2
                error('Affirmation Error: Legs which are suppose to be trace do not appear twice in this calcuation')
            end
            
            FirstLoc = min(EndLocs);
            FinalLoc = max(EndLocs);
            
            IsTraceShaded(FirstLoc:FinalLoc) = true;
            
        end
        
        ActionCyclicLegsIsTrace(Temp) = {IsTrace};
        ActionCyclicLegsIsTraceShaded(Temp) = {IsTraceShaded};
        CyclicLegsIsTraceShaded{nn} = IsTraceShaded;
        
        ActionCyclicLegsConnectedTensor(Temp) = {ConnectedTensor};
        
        ActionTensorsChargeDirectionsSides(Temp) = {ChargeSides};
        
    end
    
    
    
    ActionCyclicLegsBraidedOriginalFirstSide = ActionCyclicLegsBraidedFirstSide;
    
    for nn = 1:size(ActionTotalTensors,1)
        
        %first swap around if they are in the wrong order, this updates:
        % - ActionTotalTensors
        % - ActionTensorsLocations
        % - ActualLegsTotal
        % - ActionTensorsStructure
        % - ActionTensorsChargeDirections
        %
        % - ActionCyclicLegsBraidedFirstSide
        % - ActionCyclicLegsBraidedExternalLeg
        % - ActionCyclicLegsMadeFirst
        % - ActionCyclicLegsIsTrace
        % - ActionCyclicLegsConnectedTensor
        %
        % - TensorPathTotal
        
        
        if ActionTypeTotal(nn) ~= 1
        
            if GTPath(TensorPathTotal{nn,1},TensorPathTotal{nn,2})
                
                ActionTensorsDictionary(nn,1:2) = ActionTensorsDictionary(nn,[2,1]);
                ActionTensorsLegLocations(nn,1:2) = ActionTensorsLegLocations(nn,[2,1]);
                ActionTensorsLegEndLocations(nn,1:2) = ActionTensorsLegEndLocations(nn,[2,1]);
                ActionTensorsLegPosModified(nn,1:2) = ActionTensorsLegPosModified(nn,[2,1]);
                ActionTensorsClumptedTensors(nn,1:2) = ActionTensorsClumptedTensors(nn,[2,1]);
                
                ActionTotalTensors(nn,1:2) = ActionTotalTensors(nn,[2,1]);
                ActionTensorsLocations(nn,1:2) = ActionTensorsLocations(nn,[2,1]);
                ActualLegsTotal(nn,1:2) = ActualLegsTotal(nn,[2,1]);
                ActionTensorsStructure(nn,1:2) = ActionTensorsStructure(nn,[2,1]);
                ActionTensorsChargeDirectionsExternal(nn,1:2) = ActionTensorsChargeDirectionsExternal(nn,[2,1]);
                ActionTensorsChargeDirectionsExternalOriginal(nn,1:2) = ActionTensorsChargeDirectionsExternalOriginal(nn,[2,1]);
                ActionTensorsChargeDirections(nn,1:2) = ActionTensorsChargeDirections(nn,[2,1]);
                ActionTensorsChargeDirectionsSides(nn,1:2) = ActionTensorsChargeDirectionsSides(nn,[2,1]);
                
                ActionCyclicLegsBraidedFirstSide(nn,1:2) = ActionCyclicLegsBraidedFirstSide(nn,[2,1]);
                ActionCyclicLegsBraidedCountFlipSide(nn,1:2) = ActionCyclicLegsBraidedCountFlipSide(nn,[2,1]);
                ActionCyclicLegsBraidedOriginalFirstSide(nn,1:2) = ActionCyclicLegsBraidedOriginalFirstSide(nn,[2,1]);
                ActionCyclicLegsBraidedExternalLeg(nn,1:2) = ActionCyclicLegsBraidedExternalLeg(nn,[2,1]);
                ActionCyclicLegsIsTrace(nn,1:2) = ActionCyclicLegsIsTrace(nn,[2,1]);
                ActionCyclicLegsIsTraceShaded(nn,1:2) = ActionCyclicLegsIsTraceShaded(nn,[2,1]);
                ActionCyclicLegsConnectedTensor(nn,1:2) = ActionCyclicLegsConnectedTensor(nn,[2,1]);
                
                TensorPathTotal(nn,1:2) = TensorPathTotal(nn,[2,1]);
                DropTensor(nn,1:2) = DropTensor(nn,[2,1]);
            end
            
            %if this is a trace then it is unimportant, otherwise we can
            %work out what the updated details have to be:
            
            Temp = ActionTotalTensors == ActionTotalTensors(nn,3);
            
            NewLegPosModified = ActionTensorsLegPosModified{nn,1}|ActionTensorsLegPosModified{nn,2};
            
            NewLegLocations = ActionTensorsLegLocations{nn,1};
            NewLegLocations(NewLegPosModified) = ActionTotalTensors(nn,3);
            
            NewClumptedTensors = [ActionTensorsClumptedTensors{nn,1},ActionTensorsClumptedTensors{nn,2}];
            
            
            ActionTensorsLegLocations(Temp) = {NewLegLocations};
            ActionTensorsLegPosModified(Temp) = {NewLegPosModified};
            ActionTensorsClumptedTensors(Temp) = {NewClumptedTensors};
            
        else
            Temp = ActionTotalTensors == ActionTotalTensors(nn,3);
            
            NewLegLocations = ActionTensorsLegLocations{nn,1};
            NewLegLocations(NewLegLocations == ActionTotalTensors(nn,1)) = ActionTotalTensors(nn,3);
            ActionTensorsLegLocations(Temp) = {NewLegLocations};
            ActionTensorsLegPosModified(Temp) = ActionTensorsLegPosModified(nn,1);
            NewClumptedTensors = ActionTensorsClumptedTensors{nn,1};
            ActionTensorsClumptedTensors(Temp) = ActionTensorsClumptedTensors(nn,1);
            
        end
        
        TensorPathTotal(ActionTotalTensors == ActionTotalTensors(nn,3)) = TensorPathTotal(nn,1);
        
        %disp(['nn = ',num2str(nn)])%DebugDisp
        
        
        if ActionTypeTotal(nn) == 0 
            %Working with Kronica
            
            %first update the terms on permuting and update Legs
            % - ActualLegsTotal
            % - ActionTensorsPermute (this includes permuting at the
            % beginning)
            % - ActionTensorsStructure
            % - ActionTensorsChargeDirections
            
            %also need to work out the rotation and where we insert B in A:
            % - RotateB
            % - InsertBAfter
            %Note that we can think of this as some contraction using a
            %trivial connector.
            
            
            
            % ActionCyclicLegsBraidedFirstSide
            % ActionCyclicLegsBraidedExternalLeg
            % ActionCyclicLegsMadeFirst
            % ActionCyclicLegsIsTrace
            % ActionCyclicLegsConnectedTensor
            
            %the rotation for B will be based on where this swaps from 
            %leftwards to righwards (I have
            %assumed I've successfully built this structure correctly)
            
            %the location for inserting into A will be the same.
            Keep = ~ActionCyclicLegsExternalLeg{nn,2} & ~ActionCyclicLegsIsTraceShaded{nn,2};
            Numbers = 1:numel(ActionCyclicLegsBraidedFirstSide{nn,2});
            Numbers = Numbers(Keep);
            RotateB = find(ActionCyclicLegsBraidedFirstSide{nn,2}(Keep(1:(end-1)))~=...
                ActionCyclicLegsBraidedFirstSide{nn,2}(Keep{nn,2}(2:end)));
            if isempty(RotateB)
                RotateB = 1;
            else
                RotateB = Numbers(RotateB)+1;
            end
            
            if numel(RotateB)>1
                error('Affirmation Error: There should only be one swap from right to left in the RotateB')
            end
            
            
            %need to adjust for external legs.
            Keep = ~ActionCyclicLegsExternalLeg{nn,1} & ~ActionCyclicLegsIsTraceShaded{nn,1};
            Numbers = 1:numel(ActionCyclicLegsBraidedFirstSide{nn,2});
            Numbers = Numbers(Keep);
            InsertBAfter = find(ActionCyclicLegsBraidedFirstSide{nn,1}(Keep(1:(end-1)))~=ActionCyclicLegsBraidedFirstSide{nn,1}(Keep(2:end)));
            if isempty(InsertBAfter)
                InsertBAfter = numel(ActionCyclicLegsBraidedFirstSide{nn,1});
            else
                InsertBAfter = Numbers(InsertBAfter);
            end
            %note that this can never be 0
            
            if numel(InsertBAfter)>1
                error('Affirmation Error: There should only be one swap from right to left in InsertBAfter')
            end
            
            
            %for Kronica we just need to permute the second case and then
            %Index is just numbers 1 to NumberLegs;
            NumberLegsA = numel(ActualLegsTotal{nn,1});
            NumberLegsB =  numel(ActualLegsTotal{nn,2});
            NumberLegs = NumberLegsA+NumberLegsB;
            
            FirstKeep = 1:NumberLegsA;
            %need some sort of ordering here for SecondKeep
            SecondKeep = [RotateB:NumberLegsB, 1:(RotateB-1)];
            
            
            Permute1 = 1:NumberLegsA;
            Permute2 = [RotateB:NumberLegsB, 1:(RotateB-1)];
            Permute3 = [1:InsertBAfter,NumberLegsA+(1:NumberLegsB),(InsertBAfter+1):NumberLegsA];
            
            %define FirstKeep, SecondKeep and Index (for rearranging the
            %kept labels
            ChargeDirectionsExternalA = ActionTensorsChargeDirectionsExternal{nn,1};
            ChargeDirectionsExternalB = ActionTensorsChargeDirectionsExternal{nn,2};
            ChargeDirectionsExternalCInit = [ChargeDirectionsExternalA(FirstKeep),ChargeDirectionsExternalB(SecondKeep)];
            ChargeDirectionsExternalC = ChargeDirectionsExternalCInit(Permute3);
            
            ChargeDirectionsSidesA = ActionTensorsChargeDirectionsSides{nn,1};
            ChargeDirectionsSidesB = ActionTensorsChargeDirectionsSides{nn,2};
            ChargeDirectionsSidesCInit = [ChargeDirectionsSidesA(FirstKeep),ChargeDirectionsSidesB(SecondKeep)];
            ChargeDirectionsSidesC = ChargeDirectionsSidesCInit(Permute3);
            
            
            TempRaw = [ActionCyclicLegsBraidedExternalLeg{nn,1}(FirstKeep),ActionCyclicLegsBraidedExternalLeg{nn,2}(SecondKeep)];
            ActionCyclicLegsBraidedExternalLeg(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActionCyclicLegsIsTrace{nn,1}(FirstKeep),ActionCyclicLegsIsTrace{nn,2}(SecondKeep)];
            ActionCyclicLegsIsTrace(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActionCyclicLegsIsTraceShaded{nn,1}(FirstKeep),ActionCyclicLegsIsTraceShaded{nn,2}(SecondKeep)];
            ActionCyclicLegsIsTraceShaded(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActionCyclicLegsConnectedTensor{nn,1}(FirstKeep),ActionCyclicLegsConnectedTensor{nn,2}(SecondKeep)];
            ActionCyclicLegsConnectedTensor(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActualLegsTotal{nn,1}(FirstKeep),ActualLegsTotal{nn,2}(SecondKeep)];
            ActualLegsTotal(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            
            
            if NumberLegs < 2
                StructureTemp = [];
                ChargeDirectionsTemp = [];
            elseif NumberLegs == 2
                StructureTemp = [-1;-2];
                ChargeDirectionsTemp = [ChargeDirectionsExternalC(1);ChargeDirectionsExternalC(2)];
            else
                StructureTemp = [[NumberLegs-2;-NumberLegs],[-1;-2],[1:(NumberLegs-3);-(3:(NumberLegs-1))]];
                ChargeDirectionsTemp = [[+1;ChargeDirectionsExternalC(NumberLegs)],[ChargeDirectionsExternalC(1);ChargeDirectionsExternalC(2)],...
                    [ones([1,NumberLegs-3]);ChargeDirectionsExternalC(3:(NumberLegs-1))]];
            end
            
            
            %Now we work out the convertion information (use the
            %FirstStandard code to rotate these tensors initially and just 
            %compute the overlap).
            % - ActionTensorsConvertWords
            % - ActionTensorsCombineCharges
            if NumberLegsA*NumberLegsB == 0
                if NumberLegsA == 0 && NumberLegsB == 0
                    ActionTensorsCombineCharges{nn,1} = [];
                    ActionTensorsCombineCharges{nn,2} = [];
                    ActionTensorsCombineCharges{nn,3} = [];
                elseif NumberLegsA == 0
                    ActionTensorsCombineCharges{nn,1} = [1:NumberLegsB];
                    ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsB-2)];
                    ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsB-2)];
                else %NumberLegsB == 0
                    ActionTensorsCombineCharges{nn,1} = [1:NumberLegsA];
                    ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsA-2)];
                    ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsA-2)];
                end
            else
                if NumberLegsB>=2
                    
                    if InsertBAfter == 1
                        if 1 ~=NumberLegsA
                            ActionTensorsCombineCharges{nn,1} = [1,NumberLegsA+(1:NumberLegsB),2:NumberLegsA]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,-1,1:(NumberLegsA-2)];%Internal
                            ActionTensorsCombineCharges{nn,3} = [(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,NaN,(InsertBAfter-1+1):(NumberLegsA-2)];%Multiplicities
                        else
                            ActionTensorsCombineCharges{nn,1} = [NaN,1+(1:NumberLegsB)]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [-2,(1:(NumberLegsB-2))];%Internal
                            ActionTensorsCombineCharges{nn,3} = [NaN,(1:(NumberLegsB-2))];%Multiplicities
                        end
                    elseif InsertBAfter == NumberLegsA
                        ActionTensorsCombineCharges{nn,1} = [1:NumberLegsA,NumberLegsA+(1:NumberLegsB)]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsA-2),NaN,-(NumberLegsA+1),(NumberLegsA-2)+(1:(NumberLegsB-2))];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsA-2),NaN,NaN,(NumberLegsA-2)+(1:(NumberLegsB-2))];%Multiplicities
                    else
                        ActionTensorsCombineCharges{nn,1} = [1:InsertBAfter,NumberLegsA+(1:NumberLegsB),(InsertBAfter+1):NumberLegsA]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(InsertBAfter-1),(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,(InsertBAfter-1):(NumberLegsA-2)];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(InsertBAfter-1),(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,NaN,(InsertBAfter-1+1):(NumberLegsA-2)];%Multiplicities
                    end
                else %NumberLegsB == 1
                    
                    if InsertBAfter == 1
                        if 1 ~=NumberLegsA
                            ActionTensorsCombineCharges{nn,1} = [1,NaN,2:NumberLegsA]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [-1,1:(NumberLegsA-2)];%Internal
                            ActionTensorsCombineCharges{nn,3} = [NaN,1:(NumberLegsA-2)];%Multiplicities
                        else
                            ActionTensorsCombineCharges{nn,1} = [NaN,NaN]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [];%Internal
                            ActionTensorsCombineCharges{nn,3} = [];%Multiplicities
                        end
                    elseif InsertBAfter == NumberLegsA
                        ActionTensorsCombineCharges{nn,1} = [1:NumberLegsA,NaN]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsA-2),NaN];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsA-2),NaN];%Multiplicities
                    else
                        ActionTensorsCombineCharges{nn,1} = [1:InsertBAfter,NaN,(InsertBAfter+1):NumberLegsA]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(InsertBAfter-1),(InsertBAfter-1):(NumberLegsA-2)];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(InsertBAfter-1),NaN,(InsertBAfter-1+1):(NumberLegsA-2)];%Multiplicities
                    end
                end
            end
            
            BraidConvertA = cell([0,2]);
            
            
            
            
            
            
            
            WorkingDictionary = [ActionTensorsDictionary{nn,1};ActionTensorsDictionary{nn,2}];
            
            
            %GlobalTwists
            
            TempFirstSideB = ActionCyclicLegsBraidedFirstSide{nn,2};
            
            GlobalTwists = cell([0,2]);
            
            %now we have the twists
            for kk = 1:numel(Twists)
                Locations = sort(Twists{kk});
                
                if numel(Locations)>1
                    
                    PullOut = cell([0,2]);
                    
                    %first pull all legs inbetween out underneath and to the
                    %right
                    FirstSpot = Locations(1);
                    LastSpot = Locations(end);
                    
                    if FirstSpot+numel(Locations)-1 ~= LastSpot %then there should be something in between
                        
                        Numbers = FirstSpot:LastSpot;
                        Numbers(Locations-FirstSpot+1) = [];
                        for ll = Numbers(end:-1:1)
                            for mm = 1:sum(Numbers>ll)
                                gg = (ll+(mm-1));
                                cLabel = -gg;
                                cDir = ChargeDirectionsExternalB(gg);
                                
                                bLabel = -(gg+1);
                                bDir = ChargeDirectionsExternalB(gg+1);
                                
                                dLabel = gg;
                                dDir = +1;
                                
                                eLabel = gg-1;
                                eDir = +1;
                                MultLabel1 = eLabel;
                                MultLabel2 = dLabel;
                                
                                if gg == 1
                                    PullOut = [PullOut; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                                elseif gg == 2
                                    aLabel = -1;
                                    aDir = ChargeDirectionsExternalB(1);
                                    PullOut = [PullOut; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                                else
                                    aLabel = gg-2;
                                    aDir = +1;
                                    PullOut = [PullOut; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv2')
                                end
                            end
                        end
                    end
                    
                    
                    
                    GlobalTwists = [GlobalTwists;PullOut];
                    
                    
                    %now reverse the pull out
                    
                    PushIn = PullOut(end:-1:1,:);
                    for ll = 1:size(PushIn,1)
                        if isequal(PushIn{ll,1},'RInv')
                            PushIn{ll,1} = 'R';
                        elseif isequal(PushIn{ll,1},'BInv')
                            PushIn{ll,1} = 'B';
                        else
                            error('Affirmation Error: PullOut has an unrecognised element')
                        end
                    end
                    
                    %now perform the correct braiding
                    
                    for pp = 1:2
                        %because it always happens twice
                        
                        if TwistType(kk) == -1 %(Left)
                            MonoBraiding = [];
                            for gg = (numel(Locations)-1):-1:1
                                MonoBraiding = [MonoBraiding,(numel(Locations)-1):-1:(numel(Locations)-gg)];
                            end
                        elseif TwistType(kk) == -1 %(Right)
                            MonoBraiding = [];
                            for gg = 1:(numel(Locations)-1)
                                MonoBraiding = [MonoBraiding,1:(numel(Locations)-gg)];
                            end
                        else
                            Error('Affirmation Error: Wrong entry in twist type');
                        end
                        
                        for gg = MonoBraiding
                            cLabel = -gg;
                            cDir = ChargeDirectionsExternalB(gg);
                            
                            bLabel = -(gg+1);
                            bDir = ChargeDirectionsExternalB(gg+1);
                            
                            dLabel = gg;
                            dDir = +1;
                            
                            eLabel = gg-1;
                            eDir = +1;
                            MultLabel1 = eLabel;
                            MultLabel2 = dLabel;
                            
                            %this will only occur if it is a B type
                            if gg == 2
                                aLabel = -1;
                                aDir = ChargeDirectionsExternalB(1);
                            else
                                aLabel = gg-2;
                                aDir = +1;
                            end
                                    
                            if TwistType(kk) == -1 %(Left)
                                if gg == 1
                                    GlobalTwists = [GlobalTwists; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                                elseif gg == 2
                                    GlobalTwists = [GlobalTwists; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                                end
                            else %TwistType(kk) == -1 %(Right)
                                if gg == 1
                                    GlobalTwists = [GlobalTwists; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                                elseif gg == 2
                                    GlobalTwists = [GlobalTwists; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                                end
                            end
                        end
                        
                    end
                    
                    GlobalTwists = [GlobalTwists;PushIn];
                end
            end
            
            TempFirstSideB = TempFirstSideB(SecondKeep);
            
            TempRaw = [ActionCyclicLegsBraidedFirstSide{nn,1}(FirstKeep),TempFirstSideB];
            ActionCyclicLegsBraidedFirstSide(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            
            TempRaw = [ActionCyclicLegsBraidedOriginalFirstSide{nn,1}(FirstKeep),ActionCyclicLegsBraidedOriginalFirstSide{nn,2}(SecondKeep)];
            ActionCyclicLegsBraidedOriginalFirstSide(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            
            Beforehand = cell([0,2]); 
            if RotateB>numel(ActualLegsTotal{nn,2})
                error('Affirmation Error: RotateB is wrong');
            end
            
            %this part is twisting all that need to be pull to the left
            
            
                %this part does the rotatating
                
                
                %Fmoves to two triangles
                
            if RotateB~=1;
                if RotateB-2 > 0
                    aLabel = +(RotateB-2);
                    aDir = +1;
                elseif RotateB-2 == 0
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalB(1);
                else
                    error('Affirmation Error: Error in F-moves for rotating')
                end
                for ll = (RotateB-2):+1:(numel(ActualLegsTotal{nn,2})-3)
                    
                    if ll>(RotateB-2)
                        bLabel = +(ll);
                        bDir = +1;
                    elseif ll == (RotateB-2)
                        bLabel = -(ll+2);
                        bDir = ChargeDirectionsExternalB(ll+2);
                    else
                        error('Affirmation Error: Error in F-moves for rotating')
                    end
                    
                    cLabel = -(ll+3);
                    cDir = ChargeDirectionsExternalB(ll+3);
                    eLabel = +(ll+1);
                    eDir = +1;
                    MultLabel1 = +(ll+1);
                    
                    if ll < (numel(ActualLegsTotal{nn,2})-3)
                        dLabel = +(ll+2);
                        dDir = +1;
                        MultLabel2 = +(ll+2);
                    elseif ll == (numel(ActualLegsTotal{nn,2})-3)
                        dLabel = nan;
                        dDir = 0;
                        MultLabel2 = nan;
                    else
                        error('Affirmation Error: Error in F-moves for rotating')
                    end
                    
                    Beforehand = [Beforehand; {'F',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;...
                        aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                end
                
                LegEntryNew = ActualLegsTotal{nn,1}(InsertBAfter);
                if LegEntryNew<0
                    error('Affirmation Error: I''ve picked up a negative leg on the first tensor')
                end
                LegSideNew = ActionTensorsLegLocations{nn,1}(LegEntryNew,:) == ActionTotalTensors(nn,1);
                if sum(LegSideNew)==2
                    error('Affirmation Error: I seem to have a trace for the first tensor')
                end
                
                
                
                
                %adjust for the first triangle
                
                
                %perform the rotation at the bottom level
                
                
                
                if RotateB >2
                    bLabel = RotateB-2;
                    bDir = +1;
                else
                    bLabel = -1;
                    bDir = ChargeDirectionsExternalB(1);
                end
                
                if RotateB<numel(ActualLegsTotal{nn,2})
                    aLabel = numel(ActualLegsTotal{nn,2})-2;
                    aDir = +1;
                else
                    aLabel = -numel(ActualLegsTotal{nn,2});
                    aDir = ChargeDirectionsExternalB(numel(ActualLegsTotal{nn,2}));
                end
                
                Beforehand = [Beforehand; {'R',[aLabel,bLabel,nan,nan,nan;...
                        aDir,bDir,0,0,0]}]; %disp('R2')
                
                ChargeDirectionsExternalB = ChargeDirectionsExternalB([RotateB:end,1:(RotateB-1)]);
            
                    
                    
                    
                    
                    
                    
                    
                    
                    
                
                %prepare the swap operation
                
                PermuteExternal = [RotateB:numel(ActualLegsTotal{nn,2}),1:(RotateB-1)];
                
                if numel(ActualLegsTotal{nn,2}) > 2
                    if ~(RotateB == 2 || numel(ActualLegsTotal{nn,2}) == RotateB)
                        PermuteInternal = [(RotateB-1):(numel(ActualLegsTotal{nn,2})-2), RotateB-2, 1:(RotateB-3)];
                    else
                        %can't be both as then number ActualNumberLegs is 2
                        if RotateB == 2
                            PermuteInternal = 1:(numel(ActualLegsTotal{nn,2})-2);
                        else
                            PermuteInternal = [numel(ActualLegsTotal{nn,2})-2,1:(numel(ActualLegsTotal{nn,2})-3)];
                        end
                        
                    end
                else
                    PermuteInternal = zeros([1,0]);
                end
                
                
                [~,PermuteExternalInv] = sort(PermuteExternal);
                [~,PermuteInternalInv] = sort(PermuteInternal);
                
                Beforehand = [Beforehand; {'Swap',...
                    [-PermuteExternal,PermuteInternal;zeros(size(PermuteExternal)),zeros(size(PermuteInternal)); PermuteInternal, zeros(size(PermuteExternal));...
                        -PermuteExternalInv,PermuteInternalInv;zeros(size(PermuteExternal)),zeros(size(PermuteInternal)); PermuteInternalInv, zeros(size(PermuteExternal))]}];
                
                
                    
                    
                    
                    
                    
                    %reverse the F-moves
                    
                
                if (numel(ActualLegsTotal{nn,2})-RotateB) > 0
                    aLabel = +(numel(ActualLegsTotal{nn,2})-RotateB);
                    aDir = +1;
                elseif (numel(ActualLegsTotal{nn,2})-RotateB) == 0
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalB(1);
                else
                    error('Affirmation Error: Error in F-moves for rotating')
                end
                for ll = (numel(ActualLegsTotal{nn,2})-3):-1:(numel(ActualLegsTotal{nn,2})-RotateB)
                    
                    if ll>(numel(ActualLegsTotal{nn,2})-RotateB)
                        bLabel = +(ll);
                        bDir = +1;
                    elseif ll == (numel(ActualLegsTotal{nn,2})-RotateB)
                        bLabel = -(ll+2);
                        bDir = ChargeDirectionsExternalB(ll+2);
                    else
                        error('Affirmation Error: Error in F-moves for rotating')
                    end
                    
                    cLabel = -(ll+3);
                    cDir = ChargeDirectionsExternalB(ll+3);
                    eLabel = +(ll+1);
                    eDir = +1;
                    MultLabel1 = +(ll+1);
                    
                    if ll < (numel(ActualLegsTotal{nn,2})-3)
                        dLabel = +(ll+2);
                        dDir = +1;
                        MultLabel2 = +(ll+2);
                    elseif ll == (numel(ActualLegsTotal{nn,2})-3)
                        dLabel = nan;
                        dDir = 0;
                        MultLabel2 = nan;
                    else
                        error('Affirmation Error: Error in F-moves for rotating')
                    end
                    
                    Beforehand = [Beforehand; {'FInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;...
                        aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                end
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            BraidConvertB = cell([0,2]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This part is fixing up the reordering, going from left to
            %right, by this stage we know that the ordered contracted legs
            %are first.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%{
            BraidConvertReorderingB = cell([0,2]);
            
            
            
            
            
            WorkingConvert = cell([0,2]);
            aa = InsertBAfter;
            AA = NumberLegsA;
            BB = NumberLegsB;
            if ~(aa == AA || BB < 2 || AA<2)
                if aa ~= 1
                    aLabel = aa-1;
                    aDir = +1;
                else
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalA(1);
                end
                
                for ll = (aa+BB):-1:(aa+3)
                    bLabel = +(ll-3);
                    bDir = +1;
                    cLabel = -ll;
                    cDir = ChargeDirectionsExternalB(ll-aa);
                    dLabel = +(ll-1);
                    eLabel = +(ll-2);
                    dDir = +1;
                    eDir = +1;
                    MultLabel1 = +(ll-2);
                    MultLabel2 = +(ll-1);
                    
                    WorkingConvert = [WorkingConvert; {'FInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;...
                        aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                end
                
                ll = aa+2;
                
                bLabel = -(ll-1);
                bDir = ChargeDirectionsExternalB(ll-aa-1);
                cLabel = -ll;
                cDir = ChargeDirectionsExternalB(ll-aa);
                dLabel = +(ll-1);
                eLabel = +(ll-2);
                dDir = +1;
                eDir = +1;
                MultLabel1 = +(ll-2);
                MultLabel2 = +(ll-1);
                
                WorkingConvert = [WorkingConvert; {'FInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;...
                    aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
            end
            
            
            
            ActionTensorsConvertWords{nn,1} = BraidConvertA;
            ActionTensorsConvertWords{nn,2} = [Beforehand;BraidConvertB;BraidConvertReorderingB];
            ActionTensorsConvertWords{nn,3} = WorkingConvert;
            
            
        elseif ActionTypeTotal(nn) == 2
            %Working with Contraction
            
            %first update the terms on permuting and update Legs
            % - ActualLegsTotal
            % - ActionTensorsPermute (this includes permuting at the
            % beginning)
            % - ActionTensorsStructure
            % - ActionTensorsChargeDirections
            
            %also need to work out the rotation and where we insert B in A:
            % - RotateB
            % - InsertBAfter
            %Note that we can think of this as some contraction using a
            %trivial connector.
            
            
            
            % ActionCyclicLegsBraidedFirstSide
            % ActionCyclicLegsBraidedExternalLeg
            % ActionCyclicLegsMadeFirst
            % ActionCyclicLegsIsTrace
            % ActionCyclicLegsConnectedTensor
            
            %the rotation for B will be based on where this swaps from 
            %leftwards to righwards (I have
            %assumed I've successfully built this structure correctly)
            
            
            NumbersA = 1:size(ActualLegsTotal{nn,1},2);
            A = any(repmat(LegsContractedTotal{nn}(:),[1,size(ActualLegsTotal{nn,1},2)]) == repmat(ActualLegsTotal{nn,1},[numel(LegsContractedTotal{nn}),1]), 1);
            FirstKeep = NumbersA(~A);
            InsertBAfter = find(A,1,'last');
            NumbersA = NumbersA(A);
            if isempty(InsertBAfter)
                error('Affirmation Error: there should be at least one entry with A for contractions');
            end
            [AOrderedLogic,IndexA] = sort(A,'ascend');
            [CheckA, TrueIndexA] = sort(ActualLegsTotal{nn,1}(A),'ascend');
            TempIndexA = IndexA(AOrderedLogic);
            IndexA(AOrderedLogic) = TempIndexA(TrueIndexA);
            
            OrderedContractedA = ActualLegsTotal{nn,1}(A);
            
            
            NumbersB = 1:size(ActualLegsTotal{nn,2},2);
            B = any(repmat(LegsContractedTotal{nn}(:),[1,size(ActualLegsTotal{nn,2},2)]) == repmat(ActualLegsTotal{nn,2},[numel(LegsContractedTotal{nn}),1]), 1);
            SecondKeep = NumbersB(~B);
            RotateB = find(B,1,'first');
            SecondKeep = SecondKeep([RotateB:end,1:(RotateB-1)]);% note that there is nothing which can be thrown away before RotateB by definition.
            NumbersB = NumbersB(B);
            if isempty(RotateB)
                error('Affirmation Error: there should be at least one entry with B for contractions');
            end
            [BOrderedLogic,IndexB] = sort(B,'descend');
            [CheckB, TrueIndexB] = sort(ActualLegsTotal{nn,2}(B),'ascend');
            TempIndexB = IndexB(BOrderedLogic);
            IndexB(BOrderedLogic) = TempIndexB(TrueIndexB);
            
            OrderedContractedB = ActualLegsTotal{nn,2}(B);
            
            
            if ~isequal(CheckA, CheckB)
                error('Affirmation Failed: the contraction permutation turned out wrong')
            end

            
            NumberLegsA = numel(ActualLegsTotal{nn,1});
            NumberLegsB = numel(ActualLegsTotal{nn,2});
            NumberLegs = NumberLegsA+NumberLegsB-2*ActionNumberContractedLegs(nn);
            
            
            Permute1 = IndexA;
            Permute2 = IndexB([1:ActionNumberContractedLegs(nn),(ActionNumberContractedLegs(nn)+1):end]);
            
            %define FirstKeep, SecondKeep and Index (for rearranging the
            %kept labels
            ChargeDirectionsExternalA = ActionTensorsChargeDirectionsExternal{nn,1};
            ChargeDirectionsExternalB = ActionTensorsChargeDirectionsExternal{nn,2};
            
            
            
            %Now we work out the convertion information (use the
            %FirstStandard code to rotate these tensors initially and just 
            %compute the overlap).
            % - ActionTensorsConvertWords
            % - ActionTensorsCombineCharges
            if NumberLegsA*NumberLegsB == 0
                if NumberLegsA == 0 && NumberLegsB == 0
                    ActionTensorsCombineCharges{nn,1} = [];
                    ActionTensorsCombineCharges{nn,2} = [];
                    ActionTensorsCombineCharges{nn,3} = [];
                elseif NumberLegsA == 0
                    ActionTensorsCombineCharges{nn,1} = [1:NumberLegsB];
                    ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsB-2)];
                    ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsB-2)];
                else %NumberLegsB == 0
                    ActionTensorsCombineCharges{nn,1} = [1:NumberLegsA];
                    ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsA-2)];
                    ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsA-2)];
                end
            else
                if NumberLegsB>=2
                    
                    if InsertBAfter == 1
                        if 1 ~=NumberLegsA
                            ActionTensorsCombineCharges{nn,1} = [1,NumberLegsA+(1:NumberLegsB),2:NumberLegsA]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,-1,1:(NumberLegsA-2)];%Internal
                            ActionTensorsCombineCharges{nn,3} = [(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,NaN,(InsertBAfter-1+1):(NumberLegsA-2)];%Multiplicities
                        else
                            ActionTensorsCombineCharges{nn,1} = [NaN,1+(1:NumberLegsB)]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [-2,(1:(NumberLegsB-2))];%Internal
                            ActionTensorsCombineCharges{nn,3} = [NaN,(1:(NumberLegsB-2))];%Multiplicities
                        end
                    elseif InsertBAfter == NumberLegsA
                        ActionTensorsCombineCharges{nn,1} = [1:NumberLegsA,NumberLegsA+(1:NumberLegsB)]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsA-2),NaN,-(NumberLegsA+1),(NumberLegsA-2)+(1:(NumberLegsB-2))];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsA-2),NaN,NaN,(NumberLegsA-2)+(1:(NumberLegsB-2))];%Multiplicities
                    else
                        ActionTensorsCombineCharges{nn,1} = [1:InsertBAfter,NumberLegsA+(1:NumberLegsB),(InsertBAfter+1):NumberLegsA]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(InsertBAfter-1),(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,(InsertBAfter-1):(NumberLegsA-2)];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(InsertBAfter-1),(NumberLegsA-2)+(1:(NumberLegsB-2)),NaN,NaN,(InsertBAfter-1+1):(NumberLegsA-2)];%Multiplicities
                    end
                else %NumberLegsB == 1
                    
                    if InsertBAfter == 1
                        if 1 ~=NumberLegsA
                            ActionTensorsCombineCharges{nn,1} = [1,NaN,2:NumberLegsA]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [-1,1:(NumberLegsA-2)];%Internal
                            ActionTensorsCombineCharges{nn,3} = [NaN,1:(NumberLegsA-2)];%Multiplicities
                        else
                            ActionTensorsCombineCharges{nn,1} = [NaN,NaN]; %because we will have rotated the charges the right way by now
                            ActionTensorsCombineCharges{nn,2} = [];%Internal
                            ActionTensorsCombineCharges{nn,3} = [];%Multiplicities
                        end
                    elseif InsertBAfter == NumberLegsA
                        ActionTensorsCombineCharges{nn,1} = [1:NumberLegsA,NaN]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(NumberLegsA-2),NaN];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(NumberLegsA-2),NaN];%Multiplicities
                    else
                        ActionTensorsCombineCharges{nn,1} = [1:InsertBAfter,NaN,(InsertBAfter+1):NumberLegsA]; %because we will have rotated the charges the right way by now
                        ActionTensorsCombineCharges{nn,2} = [1:(InsertBAfter-1),(InsertBAfter-1):(NumberLegsA-2)];%Internal
                        ActionTensorsCombineCharges{nn,3} = [1:(InsertBAfter-1),NaN,(InsertBAfter-1+1):(NumberLegsA-2)];%Multiplicities
                    end
                end
            end
            
            BraidConvertA = cell([0,2]);
            StartNumber = min(NumbersA);
            EndNumber = max(NumbersA);
            for kk = (StartNumber+1):1:(EndNumber-1)
                if ~any(NumbersA(:)==kk)
                    PriorNumbers = kk+(0:-1:-(sum(NumbersA(:)<kk)-1));
                    for ll = PriorNumbers
                        if ll == 1
                            error('Affirmation Error: We should not be braiding the first one leg with anything');
                        elseif ll == 2
                            %use an R-move rather then a B-move
                            aDir = ChargeDirectionsExternalA(1);
                            bDir = ChargeDirectionsExternalA(2);
                            cDir = +1;
                            aLabel = -1;
                            bLabel = -2;
                            if numel(ActualLegsTotal{nn,1})>2
                                cLabel = +1;
                                MultLabel = +1;
                            else
                                %cLabel = NaN;
                                %MultLabel = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            ChargeDirectionsExternalA(1:2) = ChargeDirectionsExternalA([2,1]);
                            BraidConvertA = [BraidConvertA; {'R',[bLabel,aLabel,cLabel,MultLabel,MultLabel;bDir,aDir,cDir,0,0]}]; %disp('R1')
                        elseif ll == 3
                            %change the first leg of the B-move
                            aDir = ChargeDirectionsExternalA(1);
                            bDir = ChargeDirectionsExternalA(3);
                            cDir = ChargeDirectionsExternalA(2);
                            dDir = +1;
                            eDir = +1;
                            aLabel = -1;
                            bLabel = -3;
                            cLabel = -2;
                            eLabel = +1;
                            MultLabel1 = +1;
                            if numel(ActualLegsTotal{nn,1}) ~= ll
                                dLabel = +2;
                                MultLabel2 = +2;
                            else
                                %dLabel = NaN;
                                %MultLabel2 = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            ChargeDirectionsExternalA(2:3) = ChargeDirectionsExternalA([3,2]);
                            BraidConvertA = [BraidConvertA; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,...
                                      MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('B1')
                        else
                            aDir = +1;
                            bDir = ChargeDirectionsExternalA(ll);
                            cDir = ChargeDirectionsExternalA(ll-1);
                            dDir = +1;
                            eDir = +1;
                            aLabel = +(ll-3);
                            bLabel = -ll;
                            cLabel = -(ll-1);
                            eLabel = +(ll-2);
                            MultLabel1 = +(ll-2);
                            if numel(ActualLegsTotal(nn,1)) ~= ll
                                dLabel = +(ll-1);
                                MultLabel2 = +(ll-1);
                            else
                                %dLabel = NaN;
                                %MultLabel2 = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            ChargeDirectionsExternalA((ll-1):ll) = ChargeDirectionsExternalA([ll,ll-1]);
                            BraidConvertA = [BraidConvertA; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,...
                                      MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('B2')
                            
                        end
                    end
                end
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            Beforehand = cell([0,2]); 
            if RotateB>numel(ActualLegsTotal{nn,2})
                error('Affirmation Error: RotateB is wrong');
            end
            
            
            LegEntryNew = ActualLegsTotal{nn,1}(InsertBAfter);
            if LegEntryNew<0
                error('Affirmation Error: I''ve picked up a negative leg on the first tensor')
            end
            LegSideNew = ActionTensorsLegLocations{nn,1}(LegEntryNew,:) == ActionTotalTensors(nn,1);
            if sum(LegSideNew)==2
                error('Affirmation Error: I seem to have a trace for the first tensor')
            end
            
            
            
            
            
            %GlobalTwists
            FunctioningLegLocations = ActionTensorsLegLocations{nn,1};
            FunctioningLegEndLocations = ActionTensorsLegEndLocations{nn,1};
            for kk = ActionTensorsClumptedTensors{nn,2}
                FunctioningLegEndLocations(FunctioningLegLocations == kk) = ActionTensorsLegEndLocations{nn,2}(FunctioningLegLocations == kk);
                FunctioningLegLocations(FunctioningLegLocations == kk) = ActionTotalTensors(nn,2);
            end
            
            
            %FixDictionary
            WorkingDictionary = [ActionTensorsDictionary{nn,1};ActionTensorsDictionary{nn,2}];
            
            WorkingTensorPath = TensorPath;
            NewWorkingTensorPath = TensorPath;
            TensorPathA = TensorPathTotal{nn,1};
            TensorPathB = TensorPathTotal{nn,2};
            
            
            
            WantedTensors = [];
            for LegEntry = ActualLegsTotal{nn,2}
                if LegEntry>0
                    WantedTensors = [WantedTensors,FunctioningLegLocations(LegEntry,:)];
                end
            end
            
            WantedTensors(WantedTensors == ActionTotalTensors(nn,1)) = [];
            WantedTensors = unique([WantedTensors,ActionTotalTensors(nn,2)]);
            
            FixedTensors = [ActionTotalTensors(nn,1),ActionTensorsClumptedTensors{nn,1},NewOutTensorSpherical];
            for kk = unique(LegLocations(:))'
                if ~any(kk==FixedTensors)
                    if GTPath(TensorPathA,TensorPath{kk})
                        FixedTensors = [FixedTensors,kk];
                    end
                end
            end
            
            %Now add all root tensors
            for kk = 1:numel(TensorPath)
                if isequal(size(TensorPath{kk}),[2,2]) && TensorPath{kk}(2,1) == 1
                    FixedTensors = [FixedTensors,kk];
                end
            end
            
            for kk = FixedTensors
                if kk<=size(TensorPath,1)&&kk~=0
                    NewWorkingTensorPath{kk} = WorkingTensorPath{kk};
                end
            end
            
            NewFixedTensors = [FixedTensors,ActionTotalTensors(nn,2),ActionTensorsClumptedTensors{nn,2}];
            NewWantedTensors = WantedTensors;
            
            
            Counter = 1;
            while ~isempty(WantedTensors)
                
                if any(WantedTensors(Counter)==FixedTensors)
                    WantedTensors(Counter) = [];
                    if Counter>numel(WantedTensors)
                        Counter = 1;
                    end
                    continue;
                end
                
                
                Original = WantedTensors(Counter);
                if Original==ActionTotalTensors(nn,2)
                    
                    
                    %make sure we have the right FirstSide details:
                    LocationA = [];
                    LocationB = [];
                    for kk = 1:numel(LegsContractedTotal{nn})
                        LocA = find(ActualLegsTotal{nn,1}==LegsContractedTotal{nn}(kk));
                        LocB = find(ActualLegsTotal{nn,2}==LegsContractedTotal{nn}(kk));
                        if numel(LocA)~=1
                            error('Affirmation Error: We have the wrong number of entries for LocA')
                        end
                        if numel(LocB) ~=1
                            error('Affirmation Error: We have the wrong number of entries for LocB')
                        end
                        
                        LocationA = [LocationA,LocA];
                        LocationB = [LocationB,LocB];
                        
                    end
                    
                    %if the ordering of any connected Legs is disagreed on
                    %by A and B throw an error
                    
                    if any(ActionCyclicLegsBraidedFirstSide{nn,1}(LocationA) == ActionCyclicLegsBraidedFirstSide{nn,2}(LocationB))
                        [OrderingB,IndexB] = sort(LocationB);
                        Ordering = LocationA(IndexB);
                        
                        FirstLoc = find(~(ActionCyclicLegsIsTraceShaded{nn,2}|ActionCyclicLegsBraidedExternalLeg{nn,2}));
                        if isempty(FirstLoc)
                            error('Affirmation Error: FirstLoc should not be empty')
                        end
                        AnythingElse = false;
                        for kk = 1:numel(FirstLoc)
                            if ~any(FirstLoc(kk) == LocationB)
                                AnythingElse = true;
                                break;
                            end
                        end
                        %this part only updates the legs that need to be
                        %contracted, so I don't need to worry about adding
                        %to the count
                        FirstSideInfo = ActionCyclicLegsBraidedFirstSide{nn,2};
                        if isequal(Ordering,sort(Ordering,'descend'))
                            %if they are all going the same way
                            if numel(ActionTensorsClumptedTensors{nn,1}) == 1||all(ActionCyclicLegsBraidedCountFlipSide{nn,1}(LocationA)==0)
                                %then things are fine as they are because A
                                %is wrong
                            elseif numel(ActionTensorsClumptedTensors{nn,2}) == 1||all(ActionCyclicLegsBraidedCountFlipSide{nn,2}(LocationB)==0)
                                %then re-write anything that B disagree 
                                %with with what is there for A
                                FirstSideInfo(LocationB) = ~ActionCyclicLegsBraidedFirstSide{nn,1}(LocationA);
                            elseif ~AnythingElse
                                %In this case we can rotate it around to
                                %say it is last if it wasn't already
                                FirstSideInfo(LocationB) = false;
                            else
                                error('Affirmation Error: Tensors A and B are not agreeing on ordering and I can''t work it out')
                            end
                        else
                            %if they are different then:
                            First = Ordering(1)<Ordering;
                            FirstSideInfo(OrderingB(First)) = true;
                            FirstSideInfo(OrderingB(~First)) = false;
                        end
                        
                        
                        
                    else 
                        FirstSideInfo = ActionCyclicLegsBraidedFirstSide{nn,2};
                    end
                    
                    %now I can assume that FirstSide Info is good:
                    
                    
                    SafeLeg = ~(ActionCyclicLegsIsTraceShaded{nn,2});%|ActionCyclicLegsBraidedExternalLeg{nn,2});
                    if ~any(SafeLeg)
                        error('Affirmation Error: I think this was a root for B as there are no valid legs')
                    end
                    FirstLoc = find(SafeLeg,1,'first');
                    if numel(FirstLoc)~=1
                        error('Affirmation Error: LocFirst doesn''t have the right number of entries')
                    end
                    %then the source should be the first leg
                    SourceLeg = ActualLegsTotal{nn,2}(FirstLoc);
                    if SourceLeg<0
                        SourceLeg = RenameNeg(2,RenameNeg(1,:)==SourceLeg);
                    end
                    SourceTensor = FunctioningLegLocations(SourceLeg,FunctioningLegLocations(SourceLeg,:)~=Original);
                    if any(SourceLeg==RenameNeg(2,:))
                        SourceLeg = RenameNeg(1,RenameNeg(2,:)==SourceLeg);
                    end
                    
                    %work out the backwards Legs
                    LastLocs = find(SafeLeg&FirstSideInfo);
                    
                    if FirstSideInfo(FirstLoc)
                        %this means we are going backwards at this point
                        if isempty(LastLocs)
                            error('Affirmation Error: LastLocs should not be empty in this case');
                        end
                        
                        FlagCanComplete = false;
                        for BackLocation = LastLocs(end:-1:1)
                            BackLeg = ActualLegsTotal{nn,2}(BackLocation);
                            if BackLeg<0
                                BackLeg = RenameNeg(2,RenameNeg(1,:)==BackLeg);
                            end
                            BackSource = FunctioningLegLocations(BackLeg,FunctioningLegLocations(BackLeg,:)~=Original);
                            if any(BackLeg == RenameNeg(2,:))
                                BackLeg = RenameNeg(1,RenameNeg(2,:)==BackLeg);
                            end
                            
                            if BackSource == ActionTotalTensors(nn,1)
                                FlagCanComplete = true;
                                break;
                                %A can't have been fixed by this tensor
                            else
                                %work out if it is fixed by this tensor.
                                %if it is already fixed then it couldn't be
                                %fixed by this tensor
                                if any(BackSource == FixedTensors)
                                    FlagCanComplete = true;
                                    break;
                                end
                                
                                SafeLegsBack = ~CyclicLegsIsTraceShaded{BackSource};%|CyclicLegsBraidedExternalLeg{BackSource});
                                if ~any(SafeLegBack)
                                    error(['Affirmation Error: I think this was a root as there are no valid safe legs back for BackSource:, ',num2str(BackSource)])
                                end
                                FirstLocForBack = find(SafeLegBack,1,'first');
                                
                                if CyclicLegsBraidedLabelsUser{BackSource}(FirstLocForBack) == BackLeg
                                    %if this is the firstLeg then we are
                                    %building it off this leg
                                    continue;
                                else
                                    break;
                                end
                            end
                            
                        end
                        
                        if FlagCanComplete
                            %this means we are normal as the source must be
                            %first
                            if any(BackSource == FixedTensors)
                                FixedTensors = [FixedTensors,Original];
                                
                                if BackSource == ActionTotalTensors(nn,1)
                                    TensorPathB = [TensorPathA;[1,0]];
                                    Loc = find(ActualLegsTotal{nn,1}==BackLeg);
                                elseif BackSource == NewOutTensorSpherical
                                    TensorPathB = [TensorSphericalPath;[1,0]];
                                    Loc = find(TensorSphericalNegsLegs==BackLeg);
                                else
                                    TensorPathB = [WorkingTensorPath{BackSource};[1,0]];
                                    Loc = find(CyclicLegsBraidedLabelsUser{BackSource}==BackLeg);
                                end
                                if numel(Loc)~=1
                                    error('Affirmation Error: Loc should have exactly 1 value')
                                end
                                TensorPathB(end-1,2) = Loc-0.5;
                            else
                                WantedTensors = [WantedTensors, BackSource];
                            end
                        else
                            error('Affirmation Error: We should be able to Complete')
                        end
                    else
                        %this means we are normal as the source must be
                        %first
                        if any(SourceTensor == FixedTensors)
                            FixedTensors = [FixedTensors,Original];
                            
                            if SourceTensor == ActionTotalTensors(nn,1)
                                TensorPathB = [TensorPathA;[1,0]];
                                Loc = find(ActualLegsTotal{nn,1}==SourceLeg);
                            elseif SourceTensor == NewOutTensorSpherical
                                TensorPathB = [TensorSphericalPath;[1,0]];
                                Loc = find(TensorSphericalNegsLegs==SourceLeg);
                            else
                                TensorPathB = [WorkingTensorPath{SourceTensor};[1,0]];
                                Loc = find(CyclicLegsBraidedLabelsUser{SourceTensor}==SourceLeg);
                            end
                            if numel(Loc)~=1
                                error('Affirmation Error: Loc should have exactly 1 value')
                            end
                            TensorPathB(end-1,2) = Loc;
                        else
                            WantedTensors = [WantedTensors, SourceTensor];
                        end
                    end
                    
                else
                    SafeLeg = ~CyclicLegsIsTraceShaded{Original};%|CyclicLegsBraidedExternalLeg{Original});
                    if ~any(SafeLeg)
                        %if none are safe then this is an isolated part and
                        %so should be a root, therefore something is wrong
                        error(['Affirmation Error: I think this was a root as there are no valid safe legs:, ',num2str(Original)])
                        %FixedTensors = [FixedTensors, Original];
                    end
                    FirstLoc = find(SafeLeg,1,'first');
                    if numel(FirstLoc)~=1
                        error('Affirmation Error: LocFirst doesn''t have the right number of entries')
                    end
                    if CyclicLegsBraidedFirstSide{Original}(FirstLoc)
                        error('Affirmation Error: Since this is not the B tensor then this should not be first');
                    end
                    %then the source should be the first leg
                    SourceLeg = CyclicLegsBraidedLabelsUser{Original}(FirstLoc);
                    if SourceLeg<0
                        SourceLeg = RenameNeg(2,RenameNeg(1,:)==SourceLeg);
                    end
                    SourceTensor = FunctioningLegLocations(SourceLeg,FunctioningLegLocations(SourceLeg,:)~=Original);
                    if any(SourceLeg == RenameNeg(2,:))
                        SourceLeg = RenameNeg(1,RenameNeg(2,:)==SourceLeg);
                    end
                    
                    %work out what the sides truely look like:
                    FirstSide = CyclicLegsBraidedFirstSide{Original};
                    for kk = 1:numel(CyclicLegsBraidedFirstSide{Original})
                        LegEntry = CyclicLegsBraidedLabelsUser{Original}(kk);
                        
                        if LegEntry>0
                            
                            if FunctioningLegLocations(LegEntry,FunctioningLegLocations(LegEntry,:)~=Original) == ActionTotalTensors(nn,1)
                                
                                Loc = find(ActualLegsTotal{nn,1} == LegEntry);
                                if numel(Loc)~=1
                                    error('Affirmation Error: Loc doesn''t have the right number of entries')
                                end
                                FirstSide(kk) = ~ActionCyclicLegsBraidedFirstSide{nn,1}(Loc);
                                
                            elseif FunctioningLegLocations(LegEntry,FunctioningLegLocations(LegEntry,:)~=Original) == ActionTotalTensors(nn,2)
                                
                                Loc = find(ActualLegsTotal{nn,2} == LegEntry);
                                if numel(Loc)~=1
                                    error('Affirmation Error: Loc doesn''t have the right number of entries')
                                end
                                FirstSide(kk) = ~ActionCyclicLegsBraidedFirstSide{nn,2}(Loc);    
                            end
                        end
                    end
                    
                    %work out the backwards Legs
                    LastLocs = find(SafeLeg&FirstSide);
                    
                    if FirstSide(FirstLoc)
                        %this means we are going backwards at this point
                        if isempty(LastLocs)
                            error('Affirmation Error: LastLocs should not be empty in this case');
                        end
                        
                        
                        FlagCanComplete = false;
                        for BackLocation = LastLocs(end:-1:1)
                            BackLeg = CyclicLegsBraidedLabelsUser{Original}(BackLocation);
                            if BackLeg<0
                                BackLeg = RenameNeg(2,RenameNeg(1,:)==BackLeg);
                            end
                            BackSource = FunctioningLegLocations(BackLeg,FunctioningLegLocations(BackLeg,:)~=Original);
                            if any(BackLeg == RenameNeg(2,:))
                                BackLeg = RenameNeg(1,RenameNeg(2,:)==BackLeg);
                            end
                            
                            if BackSource == ActionTotalTensors(nn,1)
                                FlagCanComplete = true;
                                break;
                                %A can't have been fixed by this tensor
                            elseif BackSource == ActionTotalTensors(nn,2)
                                %if B is already fixed then it couldn't be
                                %fixed by this tensor
                                if any(BackSource == FixedTensor)
                                    FlagCanComplete = true;
                                    break;
                                end
                                
                                SafeLegsBack = ~ActionCyclicLegsIsTraceShaded{nn,2};%|ActionCyclicLegsBraidedExternalLeg{nn,2});
                                if ~any(SafeLegBack)
                                    error(['Affirmation Error: I think this was a root as there are no valid safe legs back for BackSource:, ',num2str(BackSource)])
                                end
                                FirstLocForBack = find(SafeLegBack,1,'first');
                                
                                if ActualLegsTotal{nn,2}(FirstLocForBack) == BackLeg
                                    %if this is the firstLeg then we are
                                    %building the other tensor off this leg
                                    continue;
                                else
                                    break;
                                end
                            else
                                %work out if it is fixed by this tensor.
                                %if it is already fixed then it couldn't be
                                %fixed by this tensor
                                if any(BackSource == FixedTensors)
                                    FlagCanComplete = true;
                                    break;
                                end
                                
                                SafeLegsBack = ~CyclicLegsIsTraceShaded{BackSource};%|CyclicLegsBraidedExternalLeg{BackSource});
                                if ~any(SafeLegsBack)
                                    error(['Affirmation Error: I think this was a root as there are no valid safe legs back for BackSource:, ',num2str(BackSource)])
                                end
                                FirstLocForBack = find(SafeLegsBack,1,'first');
                                
                                if CyclicLegsBraidedLabelsUser{BackSource}(FirstLocForBack) == BackLeg
                                    %if this is the firstLeg then we are
                                    %building the other tensor off this leg
                                    continue;
                                else
                                    break;
                                end
                            end
                            
                        end
                        
                        if FlagCanComplete
                            %this means we are normal as the source must be
                            %first
                            if any(BackSource == FixedTensors)
                                FixedTensors = [FixedTensors,Original];
                                
                                if BackSource == ActionTotalTensors(nn,1)
                                    WorkingTensorPath{Original} = [TensorPathA;[1,0]];
                                    Loc = find(ActualLegsTotal{nn,1}==BackLeg);
                                elseif BackSource == ActionTotalTensors(nn,2)
                                    WorkingTensorPath{Original} = [TensorPathB;[1,0]];
                                    Loc = find(ActualLegsTotal{nn,1}==BackLeg);
                                elseif BackSource == NewOutTensorSpherical
                                    WorkingTensorPath{Original} = [TensorSphericalPath;[1,0]];
                                    Loc = find(TensorSphericalNegsLegs==BackLeg);
                                else
                                    WorkingTensorPath{Original} = [WorkingTensorPath{BackSource};[1,0]];
                                    Loc = find(CyclicLegsBraidedLabelsUser{BackSource}==BackLeg);
                                end
                                if numel(Loc)~=1
                                    error('Affirmation Error: Loc should have exactly 1 value')
                                end
                                WorkingTensorPath{Original}(end-1,2) = Loc-0.5;
                            else
                                WantedTensors = [WantedTensors, BackSource];
                            end
                        else
                            error('Affirmation Error: We should be able to Complete')
                        end
                    else
                        %this means we are normal as the source must be
                        %first
                        if any(SourceTensor == FixedTensors)
                            FixedTensors = [FixedTensors,Original];
                            
                            if SourceTensor == ActionTotalTensors(nn,1)
                                WorkingTensorPath{Original} = [TensorPathA;[1,0]];
                                Loc = find(ActualLegsTotal{nn,1}==SourceLeg);
                            elseif SourceTensor == ActionTotalTensors(nn,2)
                                WorkingTensorPath{Original} = [TensorPathB;[1,0]];
                                Loc = find(ActualLegsTotal{nn,2}==SourceLeg);
                            elseif SourceTensor == NewOutTensorSpherical
                                WorkingTensorPath{Original} = [TensorSphericalPath;[1,0]];
                                Loc = find(TensorSphericalNegsLegs==SourceLeg);
                            else
                                WorkingTensorPath{Original} = [WorkingTensorPath{SourceTensor};[1,0]];
                                Loc = find(CyclicLegsBraidedLabelsUser{SourceTensor}==SourceLeg);
                            end
                            if numel(Loc)~=1
                                error('Affirmation Error: Loc should have exactly 1 value')
                            end
                            WorkingTensorPath{Original}(end-1,2) = Loc;
                        else
                            WantedTensors = [WantedTensors, SourceTensor];
                        end
                    end
                end
                
                Counter = Counter+1;
                if Counter>numel(WantedTensors)
                    Counter = 1;
                end
            end
            
            %list of locations are now LocationB and first side is
            %FirstSideInfo:
            
            
            BContractionsFirstSide = FirstSideInfo(LocationB);
            [~,IndexB] = sort(LocationB);
            BContractionsFirstSide = BContractionsFirstSide(IndexB);
            
            FlagBackwardsOnly = false;
            FlagForwardsOnly = false;
            if all(BContractionsFirstSide)
                %then we only have backwards legs
                FlagBackwardsOnly = true;
                %this means we have to twist the back to the front
                SecondKeep = 1:numel(ActualLegsTotal{nn,2});
                SecondKeep = SecondKeep([(max(LocationB)+1):end,1:max(LocationB)]);
                SecondKeep(any(repmat(SecondKeep,[numel(LocationB),1]) ==repmat(LocationB(:),[1,numel(SecondKeep)]),1)) = [];
                
                Flip = (max(LocationB)+1):numel(ActualLegsTotal{nn,2});
                %I want all these to be first
                %if any(~ActionCyclicLegsBraidedOriginalFirstSide{nn,2}(Flip))
                %    error('Error: I Haven''t coded this yet')
                %end
                FlagFlipToFront = true;
                
            elseif all(~BContractionsFirstSide)
                %then we only have forward legs
                FlagForwardsOnly = true;
                %this means we have to twist the front to the back
                SecondKeep = 1:numel(ActualLegsTotal{nn,2});
                SecondKeep = SecondKeep([(min(LocationB)):end,1:(min(LocationB)-1)]);
                SecondKeep(any(repmat(SecondKeep,[numel(LocationB),1]) ==repmat(LocationB(:),[1,numel(SecondKeep)]),1)) = [];
                
                Flip = 1:(min(LocationB)-1);
                %I want all these to be second
                %if any(ActionCyclicLegsBraidedOriginalFirstSide{nn,2}(Flip))
                %    error('Error: I Haven''t coded this yet')
                %end
                FlagFlipToFront = false;
                
            else
                %then we have a combination of the two
                %this means we have to twist the back to the front
                SecondKeep = 1:numel(ActualLegsTotal{nn,2});
                SecondKeep = SecondKeep([(max(LocationB)+1):end,1:max(LocationB)]);
                SecondKeep(any(repmat(SecondKeep,[numel(LocationB),1]) ==repmat(LocationB(:),[1,numel(SecondKeep)]),1)) = [];
                
                Flip = (max(LocationB)+1):numel(ActualLegsTotal{nn,2});
                %I want all these to be first
                if any(~ActionCyclicLegsBraidedFirstSide{nn,2}(Flip))%I don't think this has to be original
                    error('Error: I Haven''t coded this yet')
                end
                FlagFlipToFront = true;
                
            end
            SecondKeep = sort(SecondKeep);
            
            %check that all backwards are after all forwards:
            BackLocationB = LocationB(FirstSideInfo(LocationB));
            ForwardLocationB = LocationB(~FirstSideInfo(LocationB));
            
            
            if min(BackLocationB)<max(ForwardLocationB)
                error('Affirmation Error: the lines are crossing');
            end
            
            if FlagFlipToFront
                Permute3 = [1:(InsertBAfter-ActionNumberContractedLegs(nn)),...
                    NumberLegsA-ActionNumberContractedLegs(nn)+([((NumberLegsB-numel(Flip)+1):NumberLegsB)-ActionNumberContractedLegs(nn),1:(NumberLegsB-numel(Flip)-ActionNumberContractedLegs(nn))]),...
                    (InsertBAfter-ActionNumberContractedLegs(nn)+1):(NumberLegsA-ActionNumberContractedLegs(nn))];
            else
                Permute3 = [1:(InsertBAfter-ActionNumberContractedLegs(nn)),...
                    NumberLegsA-ActionNumberContractedLegs(nn)+([(numel(Flip)+1):(NumberLegsB-ActionNumberContractedLegs(nn)),1:numel(Flip)]),...
                    (InsertBAfter-ActionNumberContractedLegs(nn)+1):(NumberLegsA-ActionNumberContractedLegs(nn))];
            end
            
            ChargeDirectionsSidesA = ActionTensorsChargeDirectionsSides{nn,1};
            ChargeDirectionsSidesB = ActionTensorsChargeDirectionsSides{nn,2};
            
            ChargeDirectionsExternalCInit = [ChargeDirectionsExternalA(FirstKeep),ChargeDirectionsExternalB(SecondKeep)];
            ChargeDirectionsExternalC = ChargeDirectionsExternalCInit(Permute3);
            ChargeDirectionsSidesCInit = [ChargeDirectionsSidesA(FirstKeep),ChargeDirectionsSidesB(SecondKeep)];
            ChargeDirectionsSidesC = ChargeDirectionsSidesCInit(Permute3);
            
            TempRaw = [ActionCyclicLegsBraidedExternalLeg{nn,1}(FirstKeep),ActionCyclicLegsBraidedExternalLeg{nn,2}(SecondKeep)];
            ActionCyclicLegsBraidedExternalLeg(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActionCyclicLegsIsTrace{nn,1}(FirstKeep),ActionCyclicLegsIsTrace{nn,2}(SecondKeep)];
            ActionCyclicLegsIsTrace(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActionCyclicLegsIsTraceShaded{nn,1}(FirstKeep),ActionCyclicLegsIsTraceShaded{nn,2}(SecondKeep)];
            ActionCyclicLegsIsTraceShaded(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActionCyclicLegsConnectedTensor{nn,1}(FirstKeep),ActionCyclicLegsConnectedTensor{nn,2}(SecondKeep)];
            ActionCyclicLegsConnectedTensor(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            TempRaw = [ActualLegsTotal{nn,1}(FirstKeep),ActualLegsTotal{nn,2}(SecondKeep)];
            ActualLegsTotal(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            
            
            TempFirstSideB = ActionCyclicLegsBraidedFirstSide{nn,2};
            TempFirstSideFlipB = ActionCyclicLegsBraidedCountFlipSide{nn,2};
            if ~isempty(Flip)
                TempFirstSideB(Flip) = ~TempFirstSideB(Flip);
                TempFirstSideFlipB(Flip) = TempFirstSideFlipB(Flip)+1;
            end
            TempRaw = [ActionCyclicLegsBraidedFirstSide{nn,1}(FirstKeep),TempFirstSideB(SecondKeep)];
            ActionCyclicLegsBraidedFirstSide(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            
            TempRaw = [ActionCyclicLegsBraidedCountFlipSide{nn,1}(FirstKeep),TempFirstSideFlipB(SecondKeep)];
            ActionCyclicLegsBraidedCountFlipSide(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            
            
            TempRaw = [ActionCyclicLegsBraidedOriginalFirstSide{nn,1}(FirstKeep),ActionCyclicLegsBraidedOriginalFirstSide{nn,2}(SecondKeep)];
            ActionCyclicLegsBraidedOriginalFirstSide(ActionTotalTensors == ActionTotalTensors(nn,3)) = {TempRaw(Permute3)};
            if NumberLegs < 2
                StructureTemp = [];
                ChargeDirectionsTemp = [];
            elseif NumberLegs == 2
                StructureTemp = [-1;-2];
                ChargeDirectionsTemp = [ChargeDirectionsExternalC(1);ChargeDirectionsExternalC(2)];
            else
                StructureTemp = [[NumberLegs-2;-NumberLegs],[-1;-2],[1:(NumberLegs-3);-(3:(NumberLegs-1))]];
                ChargeDirectionsTemp = [[+1;ChargeDirectionsExternalC(NumberLegs)],[ChargeDirectionsExternalC(1);ChargeDirectionsExternalC(2)],...
                    [ones([1,NumberLegs-3]);ChargeDirectionsExternalC(3:(NumberLegs-1))]];
            end
            
            
            Counter = 1;
            %{
            while ~isempty(NewWantedTensors)
                
                if any(NewWantedTensors(Counter)==NewFixedTensors)
                    NewWantedTensors(Counter) = [];
                    if Counter>numel(NewWantedTensors)
                        Counter = 1;
                    end
                    continue;
                end
                
                
                Original = NewWantedTensors(Counter);
                
                SafeLeg = ~(CyclicLegsIsTraceShaded{Original}|CyclicLegsBraidedExternalLeg{Original});
                if ~any(SafeLeg)
                    %if none are safe then this is an isolated part and
                    %so should be a root, therefore something is wrong
                    error(['Affirmation Error: I think this was a root as there are no valid safe legs:, ',num2str(Original)])
                    %FixedTensors = [FixedTensors, Original];
                end
                FirstLoc = find(SafeLeg,1,'first');
                if numel(FirstLoc)~=1
                    error('Affirmation Error: LocFirst doesn''t have the right number of entries')
                end
                if CyclicLegsBraidedFirstSide{Original}(FirstLoc)
                    error('Affirmation Error: Since this is not the B tensor then this should not be first');
                end
                %then the source should be the first leg
                SourceLeg = CyclicLegsBraidedLabelsUser{Original}(FirstLoc);
                if SourceLeg<0
                    SourceLeg = RenameNeg(2,RenameNeg(1,:)==SourceLeg);
                end
                SourceTensor = FunctioningLegLocations(SourceLeg,FunctioningLegLocations(SourceLeg,:)~=Original);
                if any(SourceLeg==RenameNeg(2,:))
                    SourceLeg = RenameNeg(2,RenameNeg(2,:)==SourceLeg);
                end
                
                %work out what the sides truely look like:
                FirstSide = CyclicLegsBraidedFirstSide{Original};
                for kk = 1:numel(CyclicLegsBraidedFirstSide{Original})
                    LegEntry = CyclicLegsBraidedLabelsUser{Original}(kk);
                    
                    if LegEntry>0
                        
                        if any(FunctioningLegLocations(LegEntry,FunctioningLegLocations(LegEntry,:)~=Original) == ActionTotalTensors(nn,1:2))
                            
                            Loc = find(ActualLegsTotal{nn,3} == LegEntry);
                            if numel(Loc)~=1
                                error('Affirmation Error: Loc doesn''t have the right number of entries')
                            end
                            FirstSide(kk) = ~ActionCyclicLegsBraidedFirstSide{nn,3}(Loc);
                          
                        end
                    end
                end
                
                %work out the backwards Legs
                LastLocs = find(SafeLeg&FirstSide);
                
                if FirstSide(FirstLoc)
                    %this means we are going backwards at this point
                    if isempty(LastLocs)
                        error('Affirmation Error: LastLocs should not be empty in this case');
                    end
                    
                    FlagCanComplete = false;
                    for BackLocation = LastLocs(end:-1:1)
                        BackLeg = CyclicLegsBraidedLabelsUser{Original}(BackLocation);
                        BackSource = FunctioningLegLocations(BackLeg,FunctioningLegLocations(BackLeg,:)~=Original);
                        
                        if any(BackSource == ActionTotalTensors(nn,1:2))
                            FlagCanComplete = true;
                            break;
                            %A and B can't have been fixed by this tensor
                        else
                            %work out if it is fixed by this tensor.
                            %if it is already fixed then it couldn't be
                            %fixed by this tensor
                            if any(BackSource == FixedTensors)
                                FlagCanComplete = true;
                                break;
                            end
                            
                            SafeLegsBack = ~(CyclicLegsIsTraceShaded{BackSource}|CyclicLegsBraidedExternalLeg{BackSource});
                            if ~any(SafeLegBack)
                                error(['Affirmation Error: I think this was a root as there are no valid safe legs back for BackSource:, ',num2str(BackSource)])
                            end
                            FirstLocForBack = find(SafeLegBack,1,'first');
                            
                            if CyclicLegsBraidedLabelsUser{BackSource}(FirstLocForBack) == BackLeg
                                %if this is the firstLeg then we are
                                %building it off this leg
                                continue;
                            else
                                break;
                            end
                        end
                        
                    end
                    
                    if FlagCanComplete
                        %this means we are normal as the source must be
                        %first
                        if any(BackSource == NewFixedTensors)
                            NewFixedTensors = [NewFixedTensors,Original];
                           
                            if any(SourceTensor == ActionTotalTensors(nn,1:2))
                                NewWorkingTensorPath{Original} = [TensorPathA;[1,0]];
                                Loc = find(ActualLegsTotal{nn,3}==SourceLeg);
                            else
                                NewWorkingTensorPath{Original} = [NewWorkingTensorPath{SourceTensor};[1,0]];
                                Loc = find(CyclicLegsBraidedLabelsUser{SourceTensor}==SourceLeg);
                            end
                            if numel(Loc)~=1
                                error('Affirmation Error: Loc should have exactly 1 value')
                            end
                            NewWorkingTensorPath{Original}(end-1,2) = Loc-0.5;
                        else
                            NewWantedTensors = [NewWantedTensors, BackSource];
                        end
                    else
                        error('Affirmation Error: We should be able to Complete')
                    end
                else
                    %this means we are normal as the source must be
                    %first
                    if any(SourceTensor == NewFixedTensors)
                        NewFixedTensors = [FixedTensors,Original];
                        
                        if any(SourceTensor == ActionTotalTensors(nn,1:2))
                            NewWorkingTensorPath{Original} = [TensorPathA;[1,0]];
                            Loc = find(ActualLegsTotal{nn,3}==SourceLeg);
                        else
                            NewWorkingTensorPath{Original} = [NewWorkingTensorPath{SourceTensor};[1,0]];
                            Loc = find(CyclicLegsBraidedLabelsUser{SourceTensor}==SourceLeg);
                        end
                        if numel(Loc)~=1
                            error('Affirmation Error: Loc should have exactly 1 value')
                        end
                        NewWorkingTensorPath{Original}(end-1,2) = Loc;
                    else
                        NewWantedTensors = [NewWantedTensors, SourceTensor];
                    end
                end
                
                Counter = Counter+1;
                if Counter>numel(NewWantedTensors)
                    Counter = 1;
                end
            end
            %}
            
            
            
            
            
            %this part is generating a new set of LegEndLocations:
            InterestingLegs = unique([ActualLegsTotal{nn,1},ActualLegsTotal{nn,2}]);
            for LegEntry = InterestingLegs
                if LegEntry>0
                    
                    TensorNumber = FunctioningLegLocations(LegEntry,1);
                    if TensorNumber == ActionTotalTensors(nn,1)
                        LegNumber = find(ActualLegsTotal{nn,1} == LegEntry);
                        if numel(LegNumber)~=1
                            error('Affirmation Error: LegNumber has the wrong number of entries')
                        end
                        FunctioningLegEndLocations{LegEntry,1} = TensorPathA;
                        FunctioningLegEndLocations{LegEntry,1}(end,2) = LegNumber;
                        
                        if ~any(LegEntry==LegsContractedTotal{nn})
                            NewLegNumber = find(ActualLegsTotal{nn,3} == LegEntry);
                            if numel(NewLegNumber)~=1
                                error('Affirmation Error: LegNumber has the wrong number of entries')
                            end
                            %NewFunctioningLegEndLocations{LegEntry,1} = TensorPathA;
                            %NewFunctioningLegEndLocations{LegEntry,1}(end,2) = NewLegNumber;
                        end
                    elseif TensorNumber == ActionTotalTensors(nn,2)
                        LegNumber = find(ActualLegsTotal{nn,2} == LegEntry);
                        if numel(LegNumber)~=1
                            error('Affirmation Error: LegNumber has the wrong number of entries')
                        end
                        FunctioningLegEndLocations{LegEntry,1} = TensorPathB;
                        FunctioningLegEndLocations{LegEntry,1}(end,2) = LegNumber;
                        
                        if ~any(LegEntry==LegsContractedTotal{nn})
                            NewLegNumber = find(ActualLegsTotal{nn,3} == LegEntry);
                            if numel(NewLegNumber)~=1
                                error('Affirmation Error: LegNumber has the wrong number of entries')
                            end
                            %NewFunctioningLegEndLocations{LegEntry,1} = TensorPathA;
                            %NewFunctioningLegEndLocations{LegEntry,1}(end,2) = NewLegNumber;
                        end
                    else %TensorNumber is other
                        LegNumber = find(CyclicLegsBraidedLabelsUser{TensorNumber} == LegEntry);
                        if numel(LegNumber)~=1
                            error('Affirmation Error: LegNumber has the wrong number of entries')
                        end
                        NewLegNumber = LegNumber;
                        %NewFunctioningLegEndLocations{LegEntry,1} = NewWorkingTensorPath{TensorNumber};
                        %NewFunctioningLegEndLocations{LegEntry,1}(end,2) = NewLegNumber;
                        FunctioningLegEndLocations{LegEntry,1} = WorkingTensorPath{TensorNumber};
                        FunctioningLegEndLocations{LegEntry,1}(end,2) = LegNumber;
                    end
                    
                    TensorNumber = FunctioningLegLocations(LegEntry,2);
                    if TensorNumber == ActionTotalTensors(nn,1)
                        LegNumber = find(ActualLegsTotal{nn,1} == LegEntry);
                        if numel(LegNumber)~=1
                            error('Affirmation Error: LegNumber has the wrong number of entries')
                        end
                        FunctioningLegEndLocations{LegEntry,2} = TensorPathA;
                        FunctioningLegEndLocations{LegEntry,2}(end,2) = LegNumber;
                        
                        if ~any(LegEntry==LegsContractedTotal{nn})
                            NewLegNumber = find(ActualLegsTotal{nn,3} == LegEntry);
                            if numel(NewLegNumber)~=1
                                error('Affirmation Error: LegNumber has the wrong number of entries')
                            end
                            %NewFunctioningLegEndLocations{LegEntry,2} = TensorPathA;
                            %NewFunctioningLegEndLocations{LegEntry,2}(end,2) = NewLegNumber;
                        end
                    elseif TensorNumber == ActionTotalTensors(nn,2)
                        LegNumber = find(ActualLegsTotal{nn,2} == LegEntry);
                        if numel(LegNumber)~=1
                            error('Affirmation Error: LegNumber has the wrong number of entries')
                        end
                        FunctioningLegEndLocations{LegEntry,2} = TensorPathB;
                        FunctioningLegEndLocations{LegEntry,2}(end,2) = LegNumber;
                        
                        if ~any(LegEntry==LegsContractedTotal{nn})
                            NewLegNumber = find(ActualLegsTotal{nn,3} == LegEntry);
                            if numel(NewLegNumber)~=1
                                error('Affirmation Error: LegNumber has the wrong number of entries')
                            end
                            %NewFunctioningLegEndLocations{LegEntry,2} = TensorPathA;
                            %NewFunctioningLegEndLocations{LegEntry,2}(end,2) = NewLegNumber;
                        end
                    else %TensorNumber is other
                        LegNumber = find(CyclicLegsBraidedLabelsUser{TensorNumber} == LegEntry);
                        if numel(LegNumber)~=1
                            error('Affirmation Error: LegNumber has the wrong number of entries')
                        end
                        NewLegNumber = LegNumber;
                        %NewFunctioningLegEndLocations{LegEntry,2} = NewWorkingTensorPath{TensorNumber};
                        %NewFunctioningLegEndLocations{LegEntry,2}(end,2) = NewLegNumber;
                        FunctioningLegEndLocations{LegEntry,2} = WorkingTensorPath{TensorNumber};
                        FunctioningLegEndLocations{LegEntry,2}(end,2) = LegNumber;
                    end
                end
            end
            
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first step is to braid things towards the center for B
            %note that this pulls the bubbles to the outside of the
            %contracted legs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ContractedSidesBack = LocationB(FirstSideInfo(LocationB));
            ContractedSidesFront = LocationB(~FirstSideInfo(LocationB));
            
            %braid right
            MonoBraidingFront = [];
            if ~isempty(ContractedSidesFront)
                Starting = max(ContractedSidesFront);
                Ending = min(ContractedSidesFront);
                
                NonContracted = Starting:-1:Ending;
                NonContracted(any(repmat(NonContracted,[numel(ContractedSidesFront),1]) == repmat(ContractedSidesFront(:),[1,numel(NonContracted)]),1))=[];
                
                if ~isempty(NonContracted)
                    for kk = NonContracted
                        MonoBraidingFront = [MonoBraidingFront,kk-1+(1:sum(kk<ContractedSidesFront))];
                    end
                end
            end
            
            %braid left
            MonoBraidingBack = [];
            NonContractedBack = [];
            if ~isempty(ContractedSidesBack)
                Starting = min(ContractedSidesBack);
                Ending = max(ContractedSidesBack);
                
                NonContracted = Starting:1:Ending;
                NonContracted(any(repmat(NonContracted,[numel(ContractedSidesBack),1]) == repmat(ContractedSidesBack(:),[1,numel(NonContracted)]),1))=[];
                
                NonContractedBack = NonContracted;
                
                if ~isempty(NonContracted)
                    for kk = NonContracted
                        MonoBraidingBack = [MonoBraidingBack,kk-(1:sum(kk>ContractedSidesBack))];%this minus one is because we are only counting from the left
                    end
                end
            end
            
            MonoBraiding = [MonoBraidingFront,MonoBraidingBack];
            TwistType = [ones([1,numel(MonoBraidingFront)]),-ones([1,numel(MonoBraidingBack)])];
            
            MonoBraidingBubblesB = MonoBraiding;
            TwistTypeBubblesB = TwistType;
            
            BraidBubblesB = cell([0,2]);
            if ~isempty(MonoBraiding)
                for kk = 1:numel(MonoBraiding)
                    gg = MonoBraiding(kk);
                    cLabel = -gg;
                    cDir = ChargeDirectionsExternalB(gg);
                    
                    bLabel = -(gg+1);
                    bDir = ChargeDirectionsExternalB(gg+1);
                    
                    if gg~=numel(ActualLegsTotal{nn,2})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalB(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            BraidBubblesB = [BraidBubblesB; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidBubblesB = [BraidBubblesB; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            BraidBubblesB = [BraidBubblesB; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidBubblesB = [BraidBubblesB; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    ChargeDirectionsExternalB([gg,gg+1]) = ChargeDirectionsExternalB([gg+1,gg]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %second step occurs if there are any forward contractions, if
            %so then we need to flip the front Z legs to the back
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BraidFlipBForward = cell([0,2]);
            
            MonoBraiding = [];
            TwistType = [];
            
            if ~isempty(ContractedSidesFront)
                for kk = 1:(min(ContractedSidesFront)-1)
                    BraidFlipBForward = [BraidFlipBForward;...{'pUpInv',[-kk;-ChargeDirectionsExternalB(kk)]};...
                          {'RoundInv',[-kk;-ChargeDirectionsExternalB(kk)]}];
                    MonoBraiding = [MonoBraiding, 1:(numel(ActualLegsTotal{nn,2})-1)];
                    TwistType = [TwistType,ones([1,(numel(ActualLegsTotal{nn,2})-1)])];%goes under
                end
            end
            %{
            if ~isempty(Flip)
                if FlagFlipToFront
                    for kk = Flip
                        BraidFlipB = [BraidFlipB;...{'pUpInv',[-kk;-ChargeDirectionsExternalB(kk)]};...
                              {'RoundInv',[-kk;-ChargeDirectionsExternalB(kk)]}];
                        MonoBraiding = [MonoBraiding, (numel(ActualLegsTotal{nn,2})-1):-1:1];
                        TwistType = [TwistType,ones([1,(numel(ActualLegsTotal{nn,2})-1)])];
                    end
                    %now we need to pull out the rest behind all the fronts
                    
                    if ~isempty(ForwardLocationB)
                        for kk = (min(ForwardLocationB)+numel(Flip)-1):-1:1
                            MonoBraiding = [MonoBraiding, kk-1+(1:numel(ForwardLocationB))];
                            TwistType = [TwistType,ones([1,numel(ForwardLocationB)])];
                        end
                    end
                    
                else
                    %Flip was defined earlier
                    
                    for kk = Flip
                        BraidFlipB = [BraidFlipB;...{'pUp',[-kk;-ChargeDirectionsExternalB(kk)]};...
                              {'Round',[-kk;-ChargeDirectionsExternalB(kk)]}];
                        MonoBraiding = [MonoBraiding, 1:(numel(ActualLegsTotal{nn,2})-1)];
                        TwistType = [TwistType,-ones([1,(numel(ActualLegsTotal{nn,2})-1)])];
                    end
                    %now its finished;
                end
            end
            %}
            MonoBraidingFlipBForward = MonoBraiding;
            TwistTypeFlipBForward = TwistType;
            if ~isempty(MonoBraiding)
                for kk = 1:numel(MonoBraiding)
                    gg = MonoBraiding(kk);
                    cLabel = -gg;
                    cDir = ChargeDirectionsExternalB(gg);
                    
                    bLabel = -(gg+1);
                    bDir = ChargeDirectionsExternalB(gg+1);
                    
                    if gg~=numel(ActualLegsTotal{nn,2})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalB(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            BraidFlipBForward = [BraidFlipBForward; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidFlipBForward = [BraidFlipBForward; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            BraidFlipBForward = [BraidFlipBForward; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidFlipBForward = [BraidFlipBForward; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    ChargeDirectionsExternalB([gg,gg+1]) = ChargeDirectionsExternalB([gg+1,gg]);
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %third step occurs if there are any backwards contractions, if
            %so then we need to flip the back Z legs to the front
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BraidFlipBBackward = cell([0,2]);
            
            MonoBraiding = [];
            TwistType = [];
            
            if ~isempty(ContractedSidesBack)
                if ~isempty(ContractedSidesFront)
                    ExtraLegsFromFront = min(ContractedSidesFront)-1;
                else
                    ExtraLegsFromFront = 0;
                end
                for kk = numel(ActualLegsTotal{nn,2})+(0:-1:(-numel(ActualLegsTotal{nn,2})+(max(ContractedSidesBack)-ExtraLegsFromFront)+1))
                    BraidFlipBBackward = [BraidFlipBBackward;...{'pUpInv',[-kk;-ChargeDirectionsExternalB(kk)]};...
                          {'Round',[-kk;-ChargeDirectionsExternalB(kk)]}];
                    MonoBraiding = [MonoBraiding, (numel(ActualLegsTotal{nn,2})-1):-1:(numel(ContractedSidesFront)+1)];
                    TwistType = [TwistType,-ones([1,(numel(ActualLegsTotal{nn,2})-1)])];%goes under
                end
            end
            
            MonoBraidingFlipBBackward = MonoBraiding;
            TwistTypeFlipBBackward = TwistType;
            if ~isempty(MonoBraiding)
                for kk = 1:numel(MonoBraiding)
                    gg = MonoBraiding(kk);
                    cLabel = -gg;
                    cDir = ChargeDirectionsExternalB(gg);
                    
                    bLabel = -(gg+1);
                    bDir = ChargeDirectionsExternalB(gg+1);
                    
                    if gg~=numel(ActualLegsTotal{nn,2})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalB(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            BraidFlipBBackward = [BraidFlipBBackward; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidFlipBBackward = [BraidFlipBBackward; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            BraidFlipBBackward = [BraidFlipBBackward; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidFlipBBackward = [BraidFlipBBackward; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    ChargeDirectionsExternalB([gg,gg+1]) = ChargeDirectionsExternalB([gg+1,gg]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %forth step is to pull W under the backwards contracted legs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BraidLeftTwist = cell([0,2]);
            
            MonoBraiding = [];
            TwistType = [];
            
            if ~isempty(ContractedSidesBack)
                for kk = (numel(ActualLegsTotal{nn,2})-numel(ContractedSidesBack)):-1:(numel(ActualLegsTotal{nn,2})-numel(ContractedSidesBack)-numel(NonContractedBack)+1)
                    BraidLeftTwist = [BraidLeftTwist;...{'pUpInv',[-kk;-ChargeDirectionsExternalB(kk)]};...
                          {'Round',[-kk;-ChargeDirectionsExternalB(kk)]}];
                    MonoBraiding = [MonoBraiding, (numel(ActualLegsTotal{nn,2})-1-numel(ContractedSidesBack)):-1:(numel(ContractedSidesFront)+1)];
                    TwistType = [TwistType,-ones([1,(numel(ActualLegsTotal{nn,2})-1)])];%goes under
                end
            end
            
            MonoBraidingLeftTwist = MonoBraiding;
            TwistTypeLeftTwist = TwistType;
            if ~isempty(MonoBraiding)
                for kk = 1:numel(MonoBraiding)
                    gg = MonoBraiding(kk);
                    cLabel = -gg;
                    cDir = ChargeDirectionsExternalB(gg);
                    
                    bLabel = -(gg+1);
                    bDir = ChargeDirectionsExternalB(gg+1);
                    
                    if gg~=numel(ActualLegsTotal{nn,2})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalB(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            BraidLeftTwist = [BraidLeftTwist; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidLeftTwist = [BraidLeftTwist; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            BraidLeftTwist = [BraidLeftTwist; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidLeftTwist = [BraidLeftTwist; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    ChargeDirectionsExternalB([gg,gg+1]) = ChargeDirectionsExternalB([gg+1,gg]);
                end
            end
            
            
            %{
            LeftTwist = (numel(ActualLegsTotal{nn,2})-numel(ContractedSidesBack)+1):numel(ActualLegsTotal{nn,2});
            BraidLeftTwist = cell([0,2]);
            
            MonoBraiding = [];
            if ~isempty(LeftTwist)
                for kk = LeftTwist
                    BraidLeftTwist = [BraidLeftTwist;...{'pUpInv',[-kk;-ChargeDirectionsExternalB(kk)]};...
                          {'RoundInv',[-kk;-ChargeDirectionsExternalB(kk)]}];
                    MonoBraiding = [MonoBraiding, (numel(ActualLegsTotal{nn,2})-1):-1:1];
                end
            end
            TwistType = ones([1,numel(MonoBraiding)]);
            
            MonoBraidingLeftTwist = MonoBraiding;
            TwistTypeLeftTwist = TwistType;
            
            if ~isempty(MonoBraiding)
                for kk = 1:numel(MonoBraiding)
                    gg = MonoBraiding(kk);
                    cLabel = -gg;
                    cDir = ChargeDirectionsExternalB(gg);
                    
                    bLabel = -(gg+1);
                    bDir = ChargeDirectionsExternalB(gg+1);
                    
                    if gg~=numel(ActualLegsTotal{nn,2})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalB(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            BraidLeftTwist = [BraidLeftTwist; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidLeftTwist = [BraidLeftTwist; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            BraidLeftTwist = [BraidLeftTwist; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidLeftTwist = [BraidLeftTwist; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    ChargeDirectionsExternalB([gg,gg+1]) = ChargeDirectionsExternalB([gg+1,gg]);
                end
            end
            %}
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fifth step is to change all backwards contractions to forwards
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            BraidAdjustBackwards = cell([0,2]);
            
            MonoBraiding = [];
            TwistType = [];
            
            if ~isempty(ContractedSidesBack)
                for kk = (numel(ActualLegsTotal{nn,2})):-1:(numel(ActualLegsTotal{nn,2})-numel(ContractedSidesBack)+1)
                    BraidAdjustBackwards = [BraidAdjustBackwards;...{'pUpInv',[-kk;-ChargeDirectionsExternalB(kk)]};...
                          {'Round',[-kk;-ChargeDirectionsExternalB(kk)]}];
                    MonoBraiding = [MonoBraiding, (numel(ActualLegsTotal{nn,2})-1):-1:1];
                    TwistType = [TwistType,-ones([1,(numel(ActualLegsTotal{nn,2})-1)])];%goes under
                end
            end
            
            MonoBraidingAdjustBackwards = MonoBraiding;
            TwistTypeAdjustBackwards = TwistType;
            if ~isempty(MonoBraiding)
                for kk = 1:numel(MonoBraiding)
                    gg = MonoBraiding(kk);
                    cLabel = -gg;
                    cDir = ChargeDirectionsExternalB(gg);
                    
                    bLabel = -(gg+1);
                    bDir = ChargeDirectionsExternalB(gg+1);
                    
                    if gg~=numel(ActualLegsTotal{nn,2})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalB(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            BraidAdjustBackwards = [BraidAdjustBackwards; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidAdjustBackwards = [BraidAdjustBackwards; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            BraidAdjustBackwards = [BraidAdjustBackwards; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            BraidAdjustBackwards = [BraidAdjustBackwards; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    ChargeDirectionsExternalB([gg,gg+1]) = ChargeDirectionsExternalB([gg+1,gg]);
                end
            end
            
            
            
            
            
            %need to fix this term
            
            
            WorkingConvert = cell([0,2]);
            aa = InsertBAfter;
            AA = NumberLegsA;
            BB = NumberLegsB;
            if ~(aa == AA || BB < 2 || AA<2)
                if aa ~= 1
                    aLabel = aa-1;
                    aDir = +1;
                else
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalA(1);
                end
                
                for ll = (aa+BB):-1:(aa+3)
                    bLabel = +(ll-3);
                    bDir = +1;
                    cLabel = -ll;
                    cDir = ChargeDirectionsExternalB(ll-aa);
                    dLabel = +(ll-1);
                    eLabel = +(ll-2);
                    dDir = +1;
                    eDir = +1;
                    MultLabel1 = +(ll-2);
                    MultLabel2 = +(ll-1);
                    
                    WorkingConvert = [WorkingConvert; {'FInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;...
                        aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                end
                
                ll = aa+2;
                
                bLabel = -(ll-1);
                bDir = ChargeDirectionsExternalB(ll-aa-1);
                cLabel = -ll;
                cDir = ChargeDirectionsExternalB(ll-aa);
                dLabel = +(ll-1);
                eLabel = +(ll-2);
                dDir = +1;
                eDir = +1;
                MultLabel1 = +(ll-2);
                MultLabel2 = +(ll-1);
                
                WorkingConvert = [WorkingConvert; {'FInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;...
                    aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
            end
            
            
            
            
            NN = ActionNumberContractedLegs(nn);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now that we have done the F-moves we need to do RCW moves to
            %make sure everything is in a column
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for rr = 1:NN
                if aa-rr==0
                    aLabel = NaN;
                    aDir = 0;
                    bLabel = -1;
                    bDir = ChargeDirectionsExternalA(1);
                    cLabel = -1;
                    cDir = ChargeDirectionsExternalA(1);
                    MultLabel = NaN;
                else
                    if aa-rr == 1
                        aLabel = -1;
                        aDir = ChargeDirectionsExternalA(1);
                    else
                        aLabel = (aa-rr-1);
                        aDir = +1;
                    end
                    bLabel = -(aa-rr+1);
                    cLabel = (aa-rr);
                    MultLabel = (aa-rr);
                    bDir = ChargeDirectionsExternalA(aa-rr+1);
                    cDir = +1;
                end
                WorkingConvert = [WorkingConvert; {'RCW',[aLabel,bLabel,cLabel,MultLabel,MultLabel;aDir,bDir,cDir,0,0]}];
            end
            
            
            if aa>1
                MidPoint = aa-1;
            else
                MidPoint = -1;
            end
            bVec = -(aa-(0:(NN-1)));
            bPrimeVec = -(aa+(1:NN));
            
            bDirVec = ChargeDirectionsExternalA(aa-(0:(NN-1)));
            if RotateB-1+NN>NumberLegsB
                error('Affirmation Error: we shouldn''t have rotated any of the contracted legs');
            end
            bPrimeDirVec = ChargeDirectionsExternalB(mod(RotateB-1+(1:NN)-1,BB)+1);
            
            ExtraCharges = zeros([3,1]);
            ExtraDirections = zeros([2,1]);
            
            FlagTopTrivial = false;
            FlagBottomTrivial = false;
            if aa>NN+1
                cVec = +(aa-(2:(NN+1)));
                cDirVec = ones([1,NN]);
            elseif aa == NN+1
                cVec = [+(aa-(2:(NN))),-1];
                cDirVec = [ones([1,NN-1]),ChargeDirectionsExternalA(1)];
            elseif aa == NN
                %Need to fix up Extra charges for this case
                %here the top is not needed
                FlagTopTrivial = true;
                
                if aa>1
                    cVec = [+(aa-(2:(NN-1))),-1,NaN];
                    cDirVec = [ones([1,NN-2]),ChargeDirectionsExternalA(1),0];
                else
                    cVec = NaN;
                    cDirVec = 0;
                end
            else
                error('Affirmation Error: NN>aa');
            end
            
            if BB>NN
                cPrimeVec = +(aa+(0:(NN-1)));
                cPrimeDirVec = ones([1,NN]);
            elseif BB == NN
                if AA>aa
                    cPrimeVec = +(aa+(0:(NN-1)));
                    cPrimeDirVec = ones([1,NN]);
                else
                    FlagBottomTrivial = true;
                    cPrimeVec = [+(aa+(0:(NN-2))),NaN];
                    cPrimeDirVec = [ones([1,NN-1]),0];
                end
            else
                error('Affirmation Error: NN>BB');
            end
            
            if ~(FlagBottomTrivial&&FlagTopTrivial)&&(FlagBottomTrivial||FlagTopTrivial)
                %if only one of them then we might need to add additional
                %legs
                
                if FlagBottomTrivial
                    if cVec(end) == -1
                        InternalCharge = nan;
                        InternalDirection = 0;
                        ExternalCharge = nan;
                        ExternalDirection = 0;
                    elseif cVec(end) == 1
                        %note this is always before the bond so we can
                        %assume that the charges are from A.
                        InternalCharge = -1;
                        InternalDirection = ChargeDirectionsExternalA(1);
                        ExternalCharge = -2;
                        ExternalDirection = ChargeDirectionsExternalA(2);
                    else
                        InternalCharge = cVec(end)-1;
                        InternalDirection = +1;
                        ExternalCharge = bVec(end)+1;
                        ExternalDirection = ChargeDirectionsExternalA(-bVec(end)-1);
                    end
                    Multiplicity = cVec(end);
                    
                    
                    
                    ExtraCharges = [ExternalCharge; InternalCharge; Multiplicity];
                    ExtraDirections = [ExternalDirection; InternalDirection];
                else %FlagTopTrivial
                    if cPrimeVec(end) == AA+BB-2
                        InternalCharge = nan;
                        InternalDirection = 0;
                        ExternalCharge = bPrimeVec(end)-1;
                        
                        %now need to work out if this is A or B here:
                        if -(AA+BB-1) ~= bPrimeVec(end)
                            error('Affirmation Error: This doesn''t appear to be the end point')
                        end
                        if aa==AA%Loc>aa && Loc<aa+BB
                            ExternalDirection = ChargeDirectionsExternalB(mod(RotateB-2,BB)+1);%(Loc-aa);%minus 2, one for the end point and one for the mod
                        else
                            ExternalDirection = ChargeDirectionsExternalA(end);
                        end
                        
                    else
                        InternalCharge = cPrimeVec(end)+1;
                        InternalDirection = +1;
                        ExternalCharge = bPrimeVec(end)-1;
                        
                        %now need to work out if this is A or B here:
                        Loc = -bPrimeVec(end)+1;
                        if Loc>aa && Loc<aa+BB
                            ExternalDirection = ChargeDirectionsExternalB(mod(RotateB-1+(Loc-aa)-1,BB)+1);
                        else
                            if Loc ~= aa+BB
                                error('Affirmation Error: This should be at the end of B because we can''t contract past it')
                            end
                            ExternalDirection = ChargeDirectionsExternalA(Loc-BB+1);
                            %there should always be at least one more
                            %external charge as the top is trivial but the
                            %bottom isn't
                        end
                    end
                    Multiplicity = InternalCharge;
                    
                    
                    
                    ExtraCharges = [ExternalCharge; InternalCharge; Multiplicity];
                    ExtraDirections = [ExternalDirection; InternalDirection];
                end
                
                
                
            end
            
            muVec = [MidPoint,cVec(1:(end-1))];
            if muVec(end)<0
                muVec(end) = NaN;
            end
            
            muPrimeVec = cPrimeVec; %note that the NaN at the end if needed is already done
            
            
            WorkingConvert = [{'ForceSame',[[bVec;bPrimeVec];...
                                      [bDirVec;-bPrimeDirVec]]};
                                      WorkingConvert; {'Dim',[[[MidPoint;ExtraCharges;zeros([2,1])],[bVec;bPrimeVec;cVec;cPrimeVec;muVec;muPrimeVec]];...
                                      [[+1;ExtraDirections;zeros([3,1])],[bDirVec;bPrimeDirVec;cDirVec;cPrimeDirVec;zeros([2,size(muVec,2)])]]]}];
            

            
            ActionTensorsConvertWords{nn,1} = BraidConvertA;
            ActionTensorsConvertWords{nn,2} = [BraidBubblesB;BraidFlipBForward;BraidFlipBBackward;BraidLeftTwist;BraidAdjustBackwards];
            ActionTensorsConvertWords{nn,3} = WorkingConvert;
            
        else % ActionTypeTotal(nn) == 1
            %working with a trace
            
            
            %first update the terms on permuting and update Legs
            % - ActualLegsTotal
            % - ActionTensorsPermute (this includes permuting at the
            % beginning)
            % - ActionTensorsStructure
            % - ActionTensorsChargeDirections
            
            A = any(repmat(LegsContractedTotal{nn}(:),[1,size(ActualLegsTotal{nn,1},2)]) == repmat(ActualLegsTotal{nn,1},[numel(LegsContractedTotal{nn}),1]), 1);
            [AOrderedLogic,IndexA] = sort(A,'descend');
            
            
            NumbersA = 1:size(ActualLegsTotal{nn,1},2);
            NumbersA = NumbersA(A);
            [~,IndexANumbers] = sort(ActualLegsTotal{nn,1}(A),'ascend'); %pair the legs up
            NumbersA = reshape(NumbersA(IndexANumbers),[2,numel(NumbersA)/2])';
            NumbersA = NumbersA.*repmat(NumbersA(1,:)<=NumbersA(2,:),[2,1])+NumbersA([2,1],:).*repmat(NumbersA(1,:)>NumbersA(2,:),[2,1]);
            [~,IndexANumbers] = sort(NumbersA(:,1),'ascend'); %order the legs from left to right.
            NumbersA = NumbersA(IndexANumbers,:);
            
            
            
            [Check, TrueIndexA] = sort(ActualLegsTotal{nn,1}(A),'ascend');
            TempIndexA = IndexA(AOrderedLogic);
            if ~isequal(Check(1:2:end), Check(2:2:end))
                error('Affirmation Failed: the trace permutation turned out wrong')
            end
            IndexA(AOrderedLogic) = TempIndexA(TrueIndexA([1:2:end,2:2:end]));
            
            
            NumberLegs = sum(~AOrderedLogical);
            NumberLegsTotal = numel(AOrderedLogical);
            
            FirstKeep = IndexA(~AOrderedLogic);
            SecondKeep = [];
            
            Permute1 = IndexA;
            Permute2 = [];
            Permute3 = 1:NumberLegs;
            
            %define FirstKeep, SecondKeep and Index (for rearranging the
            %kept labels
            ChargeDirectionsExternalA = ActionTensorsChargeDirections{nn,1};
            ChargeDirectionsExternalC = ChargeDirectionsExternalA(FirstKeep);
            
            ChargeDirectionsSidesA = ActionTensorsChargeDirectionsSides{nn,1};
            ChargeDirectionsSidesC = ChargeDirectionsSidesA(FirstKeep);
            
            
            if NumberLegs < 2
                StructureTemp = [];
                ChargeDirectionsTemp = [];
            elseif NumberLegs == 2
                StructureTemp = [-1;-2];
                ChargeDirectionsTemp = ChargeDirectionsExternalC([1;2]);
            else
                StructureTemp = [[NumberLegs-2;-NumberLegs],[-1;-2],[1:(NumberLegs-3);-(3:(NumberLegs-1))]];
                ChargeDirectionsTemp = [[+1;ChargeDirectionsExternalC(NumberLegs)],ChargeDirectionsExternalC([1;2]),...
                    [ones([1,NumberLegs-3]);ChargeDirectionsExternalC(3:(NumberLegs-1))]];
            end
            
            ActionCyclicLegsBraidedOriginalFirstSide(ActionTotalTensors == ActionTotalTensors(nn,3)) = {ActionCyclicLegsBraidedOriginalFirstSide{nn,1}(FirstKeep)};
            
            ActualLegsTotal(ActionTotalTensors == ActionTotalTensors(nn,3)) = {[ActualLegsTotal{nn,1}(FirstKeep)]};
            
            %Now we work out the convertion information (use the
            %FirstStandard code to rotate these tensors initially and just 
            %compute the overlap).
            % - ActionTensorsConvertWords
            % - ActionTensorsCombineCharges
            
            ActionTensorsCombineCharges{nn,1} = 1:NumberLegsTotal; %because we will have rotated the charges the right way by now
            if NumberLegsTotal>2
                ActionTensorsCombineCharges{nn,2} = 1:(NumberLegsTotal-2);%Internal
                ActionTensorsCombineCharges{nn,3} = 1:(NumberLegsTotal-2);%Multiplicities
            else
                ActionTensorsCombineCharges{nn,2} = [];%Internal
                ActionTensorsCombineCharges{nn,3} = [];%Multiplicities
            end
            
            
            %now work out the words for converting
            BraidConvertA = cell([0,2]);
            StartNumber = min(NumbersA(:));
            EndNumber = max(NumbersA(:));
            for kk = (StartNumber+1):1:(EndNumber-1)
                if ~any(NumbersA(:)==kk)
                    PriorNumbers = kk+(0:-1:-(sum(NumbersA(:)<kk)-1));
                    for ll = PriorNumbers
                        if ll == 1
                            error('Affirmation Error: We should not be braiding the first one leg with anything');
                        elseif ll == 2
                            %use an R-move rather then a B-move
                            aDir = ChargeDirectionsExternalA(1);
                            bDir = ChargeDirectionsExternalA(2);
                            cDir = +1;
                            aLabel = -1;
                            bLabel = -2;
                            if numel(ActualLegsTotal(nn,1))>2
                                cLabel = +1;
                                MultLabel = +1;
                            else
                                %cLabel = NaN;
                                %MultLabel = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            ChargeDirectionsExternalA(1:2) = ChargeDirectionsExternalA([2,1]);
                            BraidConvertA = [BraidConvertA; {'R',[aLabel,bLabel,cLabel,MultLabel,MultLabel;aDir,bDir,cDir,0,0]}];
                        elseif ll == 3
                            %change the first leg of the B-move
                            aDir = ChargeDirectionsExternalA(1);
                            bDir = ChargeDirectionsExternalA(3);
                            cDir = ChargeDirectionsExternalA(2);
                            dDir = +1;
                            eDir = +1;
                            aLabel = -1;
                            bLabel = -3;
                            cLabel = -2;
                            eLabel = +1;
                            MultLabel1 = +1;
                            if numel(ActualLegsTotal{nn,1}) ~= ll
                                dLabel = +2;
                                MultLabel2 = +2;
                            else
                                %dLabel = NaN;
                                %MultLabel2 = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            ChargeDirectionsExternalA(2:3) = ChargeDirectionsExternalA([3,2]);
                            BraidConvertA = [BraidConvertA; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,...
                                      MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                        else
                            aDir = +1;
                            bDir = ChargeDirectionsExternalA(ll);
                            cDir = ChargeDirectionsExternalA(ll-1);
                            dDir = +1;
                            eDir = +1;
                            aLabel = +(ll-3);
                            bLabel = -ll;
                            cLabel = -(ll-1);
                            eLabel = +(ll-2);
                            MultLabel1 = +(ll-2);
                            if numel(ActualLegsTotal{nn,1}) ~= ll
                                dLabel = +(ll-1);
                                MultLabel2 = +(ll-1);
                            else
                                %dLabel = NaN;
                                %MultLabel2 = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            ChargeDirectionsExternalA((ll-1):ll) = ChargeDirectionsExternalA([ll,ll-1]);
                            BraidConvertA = [BraidConvertA; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,...
                                      MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                            
                        end
                    end
                end
            end
            
            %now we have removed the bubbles and just need to work out the
            %dimension contractions:
            
            
            if max(NumbersA(:)) == numel(ActionTensorsChargeDirections{nn,1})
                FlagTerminateAtRightEnd = true;
            else
                FlagTerminateAtRightEnd = false;
            end
            
            PriorList = [];
            PriorNumber = max(NumbersA(:))+1;
            PreviousUsed = 0;
            InitialPoint = max(NumbersA(:))-numel(NumbersA);%this is fixed
            
            
            for kk = 1:size(NumbersA,2)
                if NumbersA(1,kk)<PriorNumber(end)
                    PriorNumber = [PriorNumber,NumbersA(2,kk)];
                    PriorList = [PriorList,kk];
                else
                    FirstAgreement = find(PriorNumber>NumbersA(2,kk),1,'last');
                    if isempty(FirstAgreement)
                        error('Affirmation Error: there should be a non-empty First Agreement for this trace')
                    end
                    
                    UpdateDimList = PriorList(FirstAgreement:end);
                    PriorList(FirstAgreement:end) = [];
                    PriorNumber((FirstAgreement+1):end) = [];
                    
                    NumberTraced = numel(UpdateDimList);
                    
                    PriorNumber = [PriorNumber,Numbers(2,kk)];
                    PriorList = [PriorList,kk];
                    TempInitialPoint = InitialPoint+numel(PriorList);
                    MidPoint = TempInitialPoint+NumberTraced-1;
                    
                    bVec = -(TempInitialPoint+(NumberTraced:-1:1));
                    bPrimeVec = -(TempInitialPoint+NumberTraced+(1:NumberTraced));
                    
                    bDirVec = ChargeDirectionsExternalA(-bVec);
                    bPrimeDirVec = ChargeDirectionsExternalA(-bPrimeVec);
                    
                    cVec = +(TempInitialPoint+(NumberTraced:-1:1)-2);
                    cPrimeVec = +(TempInitialPoint+(1:NumberTraced)+NumberTraced-1);
                    
                    cDirVec = ones(size(bDirVec));
                    cPrimeDirVec = cDirVec;
                    
                    
                    muVec = [MidPoint,cVec(1:(end-1))];%+(TempInitialPoint+(numel(UpdateDimList):-1:1)-1);
                    muPrimeVec = cPrimeVec;%+(TempInitialPoint+(1:numel(UpdateDimList))+numel(UpdateDimList)-1);
                    
                    FlagWeirdMid = false;
                    if TempInitialPoint == 0;
                        if NumberTraced >1
                            cVec([end-1,end]) = [-1,NaN];
                            muVec(end) = NaN;
                            cDirVec([end-1,end]) = [ChargeDirectionsExternalA(1),0];
                        else
                            cVec(end) = NaN;
                            muVec(end) = NaN;
                            cDirVec(end) = 0;
                            MidPoint = -1;
                            FlagWeirdMid = ChargeDirectionsExternalA(1)==-1;
                        end
                    elseif TempInitialPoint == 1
                        cVec(end) = -1;
                        cDirVec(end) = ChargeDirectionsExternalA(1);
                    end
                    
                    PriorUsed = PriorUsed + NumberTraced;
                    ChargeDirectionsExternalA([-bVec,-bPrimeVec]) = [];
                    
                    WorkingConvert = [WorkingConvert; {'Dim',[[[MidPoint;ExtraCharges;zeros([2,1])],[bVec;bPrimeVec;cVec;cPrimeVec;muVec;muPrimeVec]];...
                                      [[+1-2*FlagWeirdMid;ExtraDirections;zeros([3,1])],[bDirVec;bDirPrimeVec;cDirVec;cDirPrimeVec;zeros([2,size(muVec,2)])]]]}];
                    
                end
            end
            
            %now contract the rest
            
            FirstAgreement = 1;
            
            UpdateDimList = PriorList(FirstAgreement:end);
            PriorList(FirstAgreement:end) = [];
            PriorNumber((FirstAgreement+1):end) = [];
            
            NumberTraced = numel(UpdateDimList);
            
            TempInitialPoint = InitialPoint;
            MidPoint = TempInitialPoint+NumberTraced-1;
            
            bVec = -(TempInitialPoint+(NumberTraced:-1:1));
            bPrimeVec = -(TempInitialPoint+NumberTraced+(1:NumberTraced));
            
            bDirVec = ChargeDirectionsExternalA(-bVec);
            bPrimeDirVec = ChargeDirectionsExternalA(-bPrimeVec);
            
            cVec = +(TempInitialPoint+(NumberTraced:-1:1)-2);
            cPrimeVec = +(TempInitialPoint+(1:NumberTraced)+NumberTraced-1);
            
            cDirVec = ones(size(bDirVec));
            cPrimeDirVec = cDirVec;
            
            
            muVec = [MidPoint,cVec(1:(end-1))];%+(TempInitialPoint+(numel(UpdateDimList):-1:1)-1);
            muPrimeVec = cPrimeVec;%+(TempInitialPoint+(1:numel(UpdateDimList))+numel(UpdateDimList)-1);
            
            FlagWeirdMid = false;
            if TempInitialPoint == 0;
                if NumberTraced >1
                    cVec([end-1,end]) = [-1,NaN];
                    muVec(end) = NaN;
                    cDirVec([end-1,end]) = [ChargeDirectionsExternalA(1),0];
                else
                    cVec(end) = NaN;
                    muVec(end) = NaN;
                    cDirVec(end) = 0;
                    MidPoint = -1;
                    FlagWeirdMid = ChargeDirectionsExternalA(1)==-1;
                end
            elseif TempInitialPoint == 1
                cVec(end) = -1;
                cDirVec(end) = ChargeDirectionsExternalA(1);
            end
            
            PriorUsed = PriorUsed + NumberTraced;
            ChargeDirectionsExternalA([-bVec,-bPrimeVec]) = [];
            
            WorkingConvert = [WorkingConvert; {'Dim',[[[MidPoint;ExtraCharges;zeros([2,1])],[bVec;bPrimeVec;cVec;cPrimeVec;muVec;muPrimeVec]];...
                                 [[+1-2*FlagWeirdMid;ExtraDirections;zeros([3,1])],[bDirVec;bPrimeDirVec;cDirVec;cPrimeDirVec;zeros([2,size(muVec,2)])]]]}];
            
            
            %need to fix this term
            ActionTensorsConvertWords{nn,1} = [];
            
            ActionTensorsConvertWords{nn,2} = [];
            
            ActionTensorsConvertWords{nn,3} = WorkingConvert;%don't need to do anything here as the dimensions are in the ActionTensorsConvertWords, this also includes 
            
            
        end
        
        
        NewLegEndLocations = ActionTensorsLegEndLocations{nn,1};
        CurrentPath = TensorPathTotal{nn,3};
        for kk = 1:numel(ActualLegsTotal{nn,3})
            FirstSide = ActionCyclicLegsBraidedFirstSide{nn,3}(kk);
            LegValue = ActualLegsTotal{nn,3}(kk);
            
            if LegValue>0
                FirstSide = ActionTensorsLegLocations{nn,1}(LegValue,:)==ActionTotalTensors(nn,1)|...
                    ActionTensorsLegLocations{nn,2}(LegValue,:)==ActionTotalTensors(nn,2);
                
                if sum(FirstSide) == 0
                    error('Affirmation Error: problem building new ActionTensorsLegEndLocations')
                elseif sum(FirstSide) == 1
                    NewLegEndLocations{LegValue,FirstSide} = CurrentPath;
                    NewLegEndLocations{LegValue,FirstSide}(end,2) = kk;
                    
                else
                    error('Affirmation Error: We shouldn''t have a contracted Leg on ActualLegsTotal{nn,3}')
                end
            end
        end
        
        ActionTensorsLegEndLocations(ActionTotalTensors == ActionTotalTensors(nn,3)) = {NewLegEndLocations};
        
        ActionTensorsStructure(ActionTotalTensors == ActionTotalTensors(nn,3)) = {StructureTemp};
        ActionTensorsChargeDirections(ActionTotalTensors == ActionTotalTensors(nn,3)) = {ChargeDirectionsTemp};
        ActionTensorsChargeDirectionsExternal(ActionTotalTensors == ActionTotalTensors(nn,3)) = {ChargeDirectionsExternalC};
        ActionTensorsChargeDirectionsExternalOriginal(ActionTotalTensors == ActionTotalTensors(nn,3)) = {ChargeDirectionsExternalC};
        ActionTensorsChargeDirectionsSides(ActionTotalTensors == ActionTotalTensors(nn,3)) = {ChargeDirectionsSidesC};
            
        ActionTensorsPermute(nn,1:3) = {Permute1,Permute2,Permute3};
        
        %working on the dictionary
        
        %first find out current position name
        %this is CurrentPath from above
        
        %next find which tensors are part of new clump
        %this is NewClumptedTensors from above
        
        
        %next find which legs are attached to the clump
        LegsOff = ActualLegsTotal{nn,3};
        
        
        
    end
    
    %Now output:
    % - C_Structure
    % - C_ChargeDirections
    
    C_Structure = ActionTensorsStructure(:,3);
    C_ChargeDirections = ActionTensorsChargeDirections(:,3);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now do the OldTrue Corrections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %need to correct:
    % - OldTrueConvert
    % - OldTruePermute
    %these three don't need to be modified, only need to make sure they
    %agree with the Convert and Permute Output
    % - OldTrueChargeSides
    % - OldTrueStructure
    % - OldTrueChargeDirections
    
    %Need to make use of:
    % - OldTrue
    
    
    BraidContractOld = cell([size(LegOrder,1),1]);
    PermuteFinal = cell([size(LegOrder,1),1]);
    %%{
    for kk = 1:size(LegOrder,1)
        
        A = ActualLegsTotal{OldTrue(kk)==ActionTotalTensors(:,3),3};
        B = LegOrder{kk,2};
        C = ActionTensorsChargeDirectionsExternal{OldTrue(kk)==ActionTotalTensors(:,3),3};
        D = ActionCyclicLegsBraidedFirstSide{OldTrue(kk)==ActionTotalTensors(:,3),3};
        
        AStore{kk} = A;
        BStore{kk} = B;
        CStore{kk} = C;
        DStore{kk} = D;
        
        OldTrueConvertRotate{kk} = cell([0,2]);
        
        if ~isempty(A)
            %if its not empty then we may need to adjust to make sure the
            %legs are in the correct sides
            
            if isempty(B)
                %if the leg order we wanted is empty and the leg order we
                %have is not empty then give an error.
                error('Affirmation Error: If the input is empty so should the output')
            end
            
            %then first perform any pDowns that are required from the
            %tensors that we are computing the environment of:
            
            if ~isempty(LegOrder{kk,1})
                LocationsEnvironments = TotalTensors(LegOrder{kk,1}(:,1))+LegOrder{kk,1}(:,2);
            else
                LocationsEnvironments = [];
            end
            
            OldTrueEnvironmentsPhase = zeros([0,2]);
            for K = LocationsEnvironments
                OldTrueEnvironmentsPhase = [OldTrueEnvironmentsPhase;ModificationsEnvironmentsPhase{K}];
            end
            
            for K = 1:size(OldTrueEnvironmentsPhase,1)
                LocLabel = find(OldTrueEnvironmentsPhase(K,1) == A);
                Dir = C(K);
                %this is because it was l bar on B while it is l on A,
                %therefore we just use the Charges on A
                OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'pDown',[-LocLabel;Dir]}];
            end
            
            %next we want to rotate to make sure Gamma is First
            
            %2 cases FlagFirst
            FlagFirst = D(1); %if this is true then the Gamma is first already, though we are going to need ajust the end points where they might not be first
            
            if FlagFirst
                Number = find(~D,1,'first');
                %first not-first leg
                if ~isempty(Number)
                    if any(D((Number+1):end))
                        error('Affirmation Error: there should be block of first and then a block of last')
                    end
                    AMod = A([Number:end,1:(Number-1)]);
                    Number = numel(D)-Number+1;
                else
                    Number = 0;
                    AMod = A;
                end
                
                MonoBraiding = [];
                for ll = 1:Number
                    K = (numel(D)-ll+1);
                    OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'Round',[-K;C(K)]}];
                    MonoBraiding = [MonoBraiding,(numel(D)-1):-1:1];
                end
                TwistType = -ones(size(MonoBraiding));%R
                
            else
                %error here
                %disp('Error Loc')
                Number = find(D,1,'first');
                %first first leg
                if ~isempty(Number)
                    if any(~D((Number+1):end))
                       error('Affirmation Error: there should be block of last and then a block of first')
                    end
                    AMod = A([Number:end,1:(Number-1)]);
                    Number = Number-1;
                else
                    Number = 0;
                    AMod = A;
                end
                
                MonoBraiding = [];
                for ll = 1:Number
                    K = ll;
                    OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'RoundInv',[-K;C(K)]}];
                    MonoBraiding = [MonoBraiding,1:(numel(D)-1)];
                end
                TwistType = ones(size(MonoBraiding));%RInv
                
            end
            
            
            
            
            if ~isempty(MonoBraiding)
                for ll = 1:numel(MonoBraiding)
                    gg = MonoBraiding(ll);
                    cLabel = -gg;
                    cDir = C(gg);
                    
                    bLabel = -(gg+1);
                    bDir = C(gg+1);
                    
                    if gg~=numel(OldTruePermute{kk})-1
                        dLabel = gg;
                        dDir = +1;
                    else
                        dLabel = nan;
                        dDir = 0;
                    end
                    
                    eLabel = gg-1;
                    eDir = +1;
                    MultLabel1 = eLabel;
                    MultLabel2 = dLabel;
                    
                    %this will only occur if it is a B type
                    if gg == 2
                        aLabel = -1;
                        aDir = C(1);
                    else
                        aLabel = gg-2;
                        aDir = +1;
                    end
                    
                    if TwistType(kk) == 1 %(Left)
                        if gg == 1
                            OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    else %TwistType(kk) == -1 %(Right)
                        if gg == 1
                            OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                        else
                            OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk}; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                        end
                    end
                    C([gg,gg+1]) = C([gg+1,gg]);
                end
            end
            
            
            %now work out required rotation:
            
            OldTrueCyclic = B(OldTrueChargeSides{kk}(1,:) == -1); 
            OldTrueCyclic = [B(OldTrueChargeSides{kk}(1,:)  == +1), OldTrueCyclic(end:-1:1)];
            %note that B and OldTrueChargeSides are both elements of LegOrder and so are uneffected by my rotation
            
            FirstLoc = find(AMod(1) == OldTrueCyclic);
            if isempty(FirstLoc)
                error('Affimation Error: didn''t find any entries for FirstLoc');
            end
            if numel(FirstLoc)  ~= 1
                error('Affimation Error: Too many entries in FirstLoc');
            end
            
            if ~isequal(OldTrueCyclic([FirstLoc:end,1:(FirstLoc-1)]), AMod)
                if FlagMustUsePlanar
                    error('Error: I can''t convert environment number %i to the form you requested', kk);
                else
                    BraidConvertOld = [BraidConvertOld,{'Fail',[]}];
                end
            end
            
            [~,IndexA] = sort(A,'ascend');
            [~,IndexB] = sort(B,'ascend');
            [~,IndexB2] = sort(IndexB,'ascend');
            
            OldTruePermute{kk}  = IndexA(IndexB2);
            
            
            %finally we want to add any Frobenius-Schur indicators when
            %type == 1
            
            if FirstLoc ~=1
                for K = (numel(D)-(FirstLoc-1)+1):numel(D)
                    OldTrueConvertRotate{kk} = [OldTrueConvertRotate{kk};{'pDown',[-K;-C(K)]}];
                end
            end
            
            
            AModStore{kk} = AMod;
            
            
        else
            if ~isempty(B)
                error('Affirmation Error: If the input is empty so should the output')
            end
        end
        Rotations(kk) = FirstLoc;
        Details{kk}{1} = OldTrueStructure{kk};
        Details{kk}{2} = OldTrueChargeDirections{kk};
        Details{kk}{3} = OldTrueChargeSides{kk}(1,:);
        Details{kk}{4} = OldTrueChargeSidesInt{kk};
        
    end
    %}
    
    %{
    for kk = 1:size(LegOrder,1)
        
        A = ActualLegsTotal{OldTrue(kk)==ActionTotalTensors(:,3),3};
        B = LegOrder{kk,2};
        C = ActionTensorsChargeDirectionsExternal{OldTrue(kk)==ActionTotalTensors(:,3),3};
        
        
        if ~isempty(A)
            
            if isempty(B)
                error('Affirmation Error: If the input is empty so should the output')
            end
            
            [~,IndexA] = sort(A,'ascend');
            [~,IndexB] = sort(B,'ascend');
            [~,IndexB2] = sort(IndexB,'ascend');
            
            OldTruePermute{kk}  = IndexA(IndexB2);
            
            
            OldTrueCyclic = B(OldTrueChargeSides{kk} == -1);
            OldTrueCyclic = [B(OldTrueChargeSides{kk}  == +1), OldTrueCyclic(end:-1:1)];
            
            FirstLoc = find(A(1) == OldTrueCyclic);
            if isempty(FirstLoc)
                error('Affimation Error: didn''t find any entries for FirstLoc');
            end
            if numel(FirstLoc)  ~= 1
                error('Affimation Error: Too many entries in FirstLoc');
            end
            
            OldTrueChargeDirectionsTemp = ActionTensorsChargeDirectionsExternal{OldTrue(kk) == ActionTotalTensors(:,3),3};
            
            %{
            if FlagAllowBraidingOldTrue
                BraidConvertOld = cell([0,2]);
                for nn = 2:(numel(OldTrueCyclic)-1)
                    nnLoc = find(A == OldTrueCyclic(mod(FirstLoc+nn-2,numel(OldTrueCyclic))+1));
                    if isempty(nnLoc)
                        error('Affimation Error: didn''t find any entries for FirstLoc');
                    end
                    if numel(nnLoc)  ~= 1
                        error('Affimation Error: Too many entries in FirstLoc');
                    end
                    
                    if nnLoc<nn
                        error('Affirmation Error: Something has gone wrong here we should be pushing unordered ones to the right so they shouldn''t pop up on the left');
                    end
                    
                    A = [A(1:(nn-1)),A(nnLoc),A(nn:(nnLoc-1)), A((nnLoc+1):end)];
                    
                    for ll = nnLoc:-1:(nn+1)
                        if ll == 1
                            error('Affirmation Error: We should not be braiding the first one leg with anything');
                        elseif ll == 2
                            %use an R-move rather then a B-move
                            aDir = OldTrueChargeDirectionsTemp(1);
                            bDir = OldTrueChargeDirectionsTemp(2);
                            cDir = +1;
                            aLabel = -1;
                            bLabel = -2;
                            if numel(A)>2
                                cLabel = +1;
                                MultLabel = +1;
                            else
                                cLabel = NaN;
                                MultLabel = NaN;
                                %error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            OldTrueChargeDirectionsTemp(1:2) = OldTrueChargeDirectionsTemp([2,1]);
                            BraidConvertOld = [BraidConvertOld; {'R',[aLabel,bLabel,cLabel,MultLabel,MultLabel;aDir,bDir,cDir,0,0]}];
                        elseif ll == 3
                            %change the first leg of the B-move
                            aDir = OldTrueChargeDirectionsTemp(1);
                            bDir = OldTrueChargeDirectionsTemp(3);
                            cDir = OldTrueChargeDirectionsTemp(2);
                            dDir = +1;
                            eDir = +1;
                            aLabel = -1;
                            bLabel = -3;
                            cLabel = -2;
                            eLabel = +1;
                            MultLabel1 = +1;
                            if numel(A) ~= ll
                                dLabel = +2;
                                MultLabel2 = +2;
                            else
                                dLabel = NaN;
                                MultLabel2 = NaN;
                                %error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            OldTrueChargeDirectionsTemp(2:3) = OldTrueChargeDirectionsTemp([3,2]);
                            BraidConvertOld = [BraidConvertOld; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,...
                                      MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                        else
                            aDir = +1;
                            bDir = OldTrueChargeDirectionsTemp(ll);
                            cDir = OldTrueChargeDirectionsTemp(ll-1);
                            dDir = +1;
                            eDir = +1;
                            aLabel = +(ll-3);
                            bLabel = -ll;
                            cLabel = -(ll-1);
                            eLabel = +(ll-2);
                            MultLabel1 = +(ll-2);
                            if numel(A) ~= ll
                                dLabel = +(ll-1);
                                MultLabel2 = +(ll-1);
                            else
                                %dLabel = NaN;
                                %MultLabel2 = NaN;
                                error('Affirmation Error: We shouldn''t have to braid the last element')
                            end
                            OldTrueChargeDirectionsTemp((ll-1):ll) = OldTrueChargeDirectionsTemp([ll,ll-1]);
                            BraidConvertOld = [BraidConvertOld; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,...
                                      MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}];
                            
                        end
                    end
                end
            else
                BraidConvertOld = cell([0,2]);
                if ~isequal(OldTrueCyclic([FirstLoc:end,1:(FirstLoc-1)]), A)
                    if FlagMustUsePlanar
                        error('Error: I can''t convert environment number %i to the form you requested', kk);
                    else
                        BraidConvertOld = [BraidConvertOld,{'Fail',[]}];
                    end
                end
            end
            %}
            
            BraidConvertOld = cell([0,2]);
            if ~isequal(OldTrueCyclic([FirstLoc:end,1:(FirstLoc-1)]), A)
                if FlagMustUsePlanar
                    error('Error: I can''t convert environment number %i to the form you requested', kk);
                else
                    BraidConvertOld = [BraidConvertOld,{'Fail',[]}];
                end
            end
            
            OldTrueConvertRotate{kk} = BraidConvertOld;
        else
            if ~isempty(B)
                error('Affirmation Error: If the input is empty so should the output')
            end
        end
        Rotations(kk) = FirstLoc;
        Details{kk}{1} = OldTrueStructure{kk};
        Details{kk}{2} = OldTrueChargeDirections{kk};
        Details{kk}{3} = OldTrueChargeSides{kk};
        Details{kk}{4} = OldTrueChargeSidesInt{kk};
        
        %now adjust to account for other things:
        
        %account for actual Frobenius Shur indicators I have to put in
        %because they are Type == 1
        if FirstLoc ~=1
            for K = 1:min(FirstLoc-1,sum(OldTrueChargeSides{kk}==+1))
                OldTrueConvertRotate{kk} = [{'pDown',[K;-C(K)]}; OldTrueConvertRotate{kk}];
            end
        end
        
        %now account for all environments
        
        %the Locations are TotalTensors(LegOrder{kk,1}(:,1))+LegOrder{kk,1}(:,2):
        
        if ~isempty(LegOrder{kk,1})
            LocationsEnvironments = TotalTensors(LegOrder{kk,1}(:,1))+LegOrder{kk,1}(:,2);
        else
            LocationsEnvironments = [];
        end
        
        OldTrueEnvironmentsPhase = zeros([0,2]);
        for K = LocationsEnvironments
            OldTrueEnvironmentsPhase = [OldTrueEnvironmentsPhase;ModificationsEnvironmentsPhase{K}];
        end
        
        for K = 1:size(OldTrueEnvironmentsPhase,1)
            LocLabel = find(OldTrueEnvironmentsPhase(K,1) == B);
            Dir = C(K);
            %this is because it was l bar on B while it is l on A,
            %therefore we just use the Charges on A
            OldTrueConvertRotate{kk} = [{'pDown',[-LocLabel;Dir]}; OldTrueConvertRotate{kk}];
        end
        
    end
    %}
    
    OldTrueConvertWordsBack = SymTensor.FromStandard(Details,Rotations, false);
    
    OldTrueConvertWords = SymTensor.FromStandard(Details,Rotations, true, false);
    
    for kk = 1:numel(OldTrueConvertRotate)
        OldTrueConvertWords{kk} = [OldTrueConvertRotate{kk};OldTrueConvertWords{kk}];
    end
    
    % now we want to update:
    
    ActionTensorsPermuteOld = ActionTensorsPermute;
    for kk = 1:length(RotationModifiersInitial)
        [I,J] = find(ActionTotalTensors == kk);
        for ll = 1:numel(I)
            %CyclicLegsBraidedTensorDetails{kk}
            
            ActionTensorsConvertWords{I(ll),J(ll)} = [RotationModifiersInitial{kk};ActionTensorsConvertWords{I(ll),J(ll)}];
            %ActionTensorsCombineCharges{I(ll),J(ll)} = ActionTensorsCombineCharges{I(ll),J(ll)};
            %This has already been done after permuting and ConvertWords,
            %so I don't need to care about it
            %ActionTensorsPermute{I(ll),J(ll)} = mod(ActionTensorsPermute{I(ll),J(ll)}-RotateTensorInitial(kk),numel(ActionTensorsPermute{I(ll),J(ll)}))+1;
            Temp = CyclicLegsBraidedTensorDetails{kk}([RotateTensorInitial(kk):end, 1:(RotateTensorInitial(kk)-1)]);
            ActionTensorsPermute{I(ll),J(ll)} = Temp(ActionTensorsPermute{I(ll),J(ll)});
        end
    end
    
end


if FlagLoadSaveBroad
    [ActionTensorsLocations, ActionTypeTotal, ActionTensorsPermute, ActionNumberContractedLegs,...
                 OldTrueLocations, OldTruePermute, OldTrueChargeSides, OldTrueChargeSidesInt,...
                 OldTrueStructure,OldTrueChargeDirections,ActionTotalTensors,DropTensor,OldTrue,OldTrueBraiding,OldTrueBraidingDirection] = SymmetryUsed.LoadCurrentData(FlagNumber);
end


if FlagNewMidDetails
    
    if FlagLoadSaveBroad
        [TotalTensors, ActualLegsTotal, LegOrder,ActionTensorsStructure,ActionTensorsChargeDirections] = SymmetryUsed.LoadDataForMidWithoutBroad(FlagNumber);
    end
    
    if FlagLoadSavePlanar
        [OldTrueConvertWords, ActionTensorsConvertWords, ActionTensorsCombineCharges,ActionTensorsChargeDirectionsExternalOriginal] = SymmetryUsed.LoadDataForMidWithoutPlanar(FlagNumber);
    end
    
    %if ~isequal(SymmetryUsed.getSymName,'none')
        ActionTensorsChargesExternal = cell(size(ActionTotalTensors));
        ActionTensorsChargesInternal = cell(size(ActionTotalTensors));
        ActionTensorsMultiplicitiesInternal = cell(size(ActionTotalTensors));
        
        for nn = 1:(numel(TotalTensors)-1)
            Temp = ActionTotalTensors>TotalTensors(nn)&ActionTotalTensors<=TotalTensors(nn+1);
            ActionTensorsChargesExternal(Temp) = repmat(ExternalChargeTensorDetails(nn), [sum(sum(Temp)),1]);
            ActionTensorsChargesInternal(Temp) = repmat(InternalChargeTensorDetails(nn), [sum(sum(Temp)),1]);
            ActionTensorsMultiplicitiesInternal(Temp) = repmat(MultiplicitiesInternalTensorDetails(nn), [sum(sum(Temp)),1]);
        end
    %end
    
    
    C_MultiplicitiesInternal = cell([size(ActionTotalTensors,1),1]);
    C_ChargeLabelsInternal = cell([size(ActionTotalTensors,1),1]);
    C_ChargeLabelsExternal = cell([size(ActionTotalTensors,1),1]);
    
    C_List = cell([size(ActionTotalTensors,1),1]);
    ActionMaxNumTensorParts = zeros([size(ActionTotalTensors,1),1]);
    ActionChargesList = cell([size(ActionTotalTensors,1),1]);
    
    
    DimensionsExternal = SymmetryUsed.getDimExternal;
    Dimensions = SymmetryUsed.getDim;
    
    TrivialIrrep = SymmetryUsed.getTrivialIrrep;
    InverseIrrep = SymmetryUsed.getInverseIrrep;
    
    if ~SymmetryUsed.IsNonAbelian && ~SymmetryUsed.IsBraiding && FlagSpeedUpAbelian
        for nn = 1:size(ActionTotalTensors,1)
            
            if ActionTypeTotal(nn) == 1 %trace
                ChargeExtTempA = ActionTensorsChargesExternal{nn,1}(:,ActionTensorsPermute{nn,1});
                ChargeExtTempA = ChargeExtTempA(:,(2*ActionNumberContractedLegs(nn)+1):end);
                [ChargeExtTempA,~,ListA] = unique(ChargeExtTempA, 'rows');
                
                ChargeExtTempC = ChargeExtTempA;
                
                if ~isempty(ActionTensorsPermute{nn,3})
                    ChargeExtTempC = ChargeExtTempC(:,ActionTensorsPermute{nn,3});
                end
                
                ActionChargesList{nn} = [ListA, zeros(size(ListA)),ListA];
                
            elseif ActionTypeTotal(nn) == 2 %contraction
                ChargeDirExtTempA = ActionTensorsChargeDirectionsExternalOriginal{nn,1}(:,ActionTensorsPermute{nn,1});
                ChargeExtTempA = ActionTensorsChargesExternal{nn,1}(:,ActionTensorsPermute{nn,1});
                ChargeExtTempA2 = ChargeExtTempA(:,1:(end-ActionNumberContractedLegs(nn)));
                [ChargeExtTempA2,~,ListA] = unique(ChargeExtTempA2, 'rows');
                
                
                ChargeDirExtTempB = ActionTensorsChargeDirectionsExternalOriginal{nn,2}(:,ActionTensorsPermute{nn,2});
                ChargeExtTempB = ActionTensorsChargesExternal{nn,2}(:,ActionTensorsPermute{nn,2});
                ChargeExtTempB2 = ChargeExtTempB(:,(ActionNumberContractedLegs(nn)+1):end);
                [ChargeExtTempB2,~,ListB] = unique(ChargeExtTempB2, 'rows');
                
                ChargeExtTempA3 = ChargeExtTempA(:,(end-ActionNumberContractedLegs(nn)+1):end);
                ChargeDirExtTempA3  = ChargeDirExtTempA(:,(end-ActionNumberContractedLegs(nn)+1):end);
                ChargeExtTempB3 = ChargeExtTempB(:,1:ActionNumberContractedLegs(nn));
                ChargeDirExtTempB3 = ChargeDirExtTempB(:,1:ActionNumberContractedLegs(nn));
                
                
                
                InverseIrrepUsed = SymmetryUsed.InverseIrrep;
                
                ChargeExtTempB3(repmat(ChargeDirExtTempB3==ChargeDirExtTempA3,[size(ChargeExtTempB3,1),1])) = ....
                    InverseIrrepUsed(ChargeExtTempB3(repmat(ChargeDirExtTempB3==ChargeDirExtTempA3,[size(ChargeExtTempB3,1),1])));
                
                
                
                ChargeExtTempC = [reshape(repmat(reshape(ChargeExtTempA2,[1,size(ChargeExtTempA2,1),size(ChargeExtTempA2,2)]),[size(ChargeExtTempB2,1),1]),...
                                          [size(ChargeExtTempA2,1)*size(ChargeExtTempB2,1),size(ChargeExtTempA2,2)]), repmat(ChargeExtTempB2,[size(ChargeExtTempA2,1),1])];
                                      
                if ~isempty(ActionTensorsPermute{nn,3})
                    ChargeExtTempC = ChargeExtTempC(:,ActionTensorsPermute{nn,3});
                end                      
                
                ActionChargesListTemp = [reshape(repmat(ListA',[size(ListB,1),1]),[size(ListA,1)*size(ListB,1),1]), repmat(ListB,[size(ListA,1),1])];
                ActionChargesList{nn} = [reshape(repmat((1:size(ListA,1)),[size(ListB,1),1]),[size(ListA,1)*size(ListB,1),1]), repmat((1:size(ListB,1))',[size(ListA,1),1])];
                ActionChargesList{nn} = [ActionChargesList{nn}, (ActionChargesListTemp(:,1)-1)*size(ChargeExtTempB2,1)+ActionChargesListTemp(:,2)];
                
                
                
                AData = ChargeExtTempA3(ActionChargesList{nn}(:,1),:);
                BData = ChargeExtTempB3(ActionChargesList{nn}(:,2),:);
                
                Keep = all(AData==BData,2);
                ActionChargesList{nn}(~Keep,:) = [];
                [Output,~,Ordering] = unique(ActionChargesList{nn}(:,3));
                ActionChargesList{nn}(:,3) = Ordering;
                ChargeExtTempC = ChargeExtTempC(Output,:);
                
            else %ifActionTypeTotal(nn) == 0 %Kronica
                ChargeExtTempA = ExternalChargeTensorDetails{nn,1}(:,ActionTensorsPermute{nn,1});
                ChargeExtTempB = ExternalChargeTensorDetails{nn,2}(:,ActionTensorsPermute{nn,2});
                ListA = (1:size(ChargeExtTempA,1))';
                ListB = (1:size(ChargeExtTempB,1))';
                
                ChargeExtTempC = [reshape(repmat(reshape(ChargeExtTempA,[1,size(ChargeExtTempA,1),size(ChargeExtTempA,2)]),[size(ChargeExtTempB,1),1]),[size(ChargeExtTempA,1)*size(ChargeExtTempB,1),size(ChargeExtTempA,2)]),...
                    repmat(ChargeExtTempB,[size(ChargeExtTempA,1),1])];
                
                if ~isempty(ActionTensorsPermute{nn,3})
                    ChargeExtTempC = ChargeExtTempC(:,ActionTensorsPermute{nn,3});
                end
                
                ActionChargesList{nn} = [reshape(repmat(ListA',[size(ListB,1),1]),[size(ListA,1)*size(ListB,1),1]), repmat(ListB,[size(ListA,1),1]), (1:(size(ListA,1)*size(ListB,1)))'];
                
            end
            
            if size(ChargeExtTempC,2) == 0
                %then there is only one output (the numeric output)
                ChargeIntTempC = zeros([1,0]);
                MultiplictiesIntTempC = zeros([1,0]);
                if size(ChargeExtTempC,1) ~=1
                    error('Affirmation Failed: ChargeExtTempC is not one entry')
                end
                
            elseif size(ChargeExtTempC,2) == 1
                ChargeIntTempC = zeros([size(ChargeExtTempC,1),0]);
                MultiplictiesIntTempC = zeros([size(ChargeExtTempC,1),0]);
                
                
                if (any(ChargeExtTempC(:,1) ~= InverseIrrep(ChargeExtTempC(:,2)))&&prod(C_ChargeDirections{nn})==1) || ...
                        (any(ChargeExtTempC(:,1) ~= ChargeExtTempC(:,2))&&prod(C_ChargeDirections{nn})==-1)
                    error('Affirmation Failed: ChargeExtTempC doesn''t have Inverse Irreps')
                end
                
            else
                ChargeDirExtTempC = ActionTensorsChargeDirectionsExternalOriginal{nn,3};%(:,ActionTensorsPermute{nn,3});
                
                ChargeExtTempCAdjust = ChargeExtTempC;
                
                ChargeExtTempCAdjust(repmat(ChargeDirExtTempC==-1,[size(ChargeExtTempCAdjust,1),1])) = ....
                    InverseIrrepUsed(ChargeExtTempCAdjust(repmat(ChargeDirExtTempC==-1,[size(ChargeExtTempCAdjust,1),1])));
                
                
                [ChargeIntTempC, MultList, MultiplicitiesIntTempC] = SymmetryUsed.FuseChargeList(ChargeExtTempCAdjust(:,1), ChargeExtTempCAdjust(:,2));
                
                NumberList = 1:size(ChargeExtTempC,1);
                DropChargesList = NumberList(~any(repmat(NumberList,[size(MultList,1),1])==repmat(MultList,[1,size(NumberList,2)]),1));
                
                if ~isempty(DropChargesList)
                    ActionChargesList{nn}(any(repmat(ActionChargesList{nn}(:,3),[1,numel(DropChargesList)])==repmat(DropChargesList,[size(ActionChargesList{nn},1),1]),2),:) = [];
                    %ActionChargesList{nn}(:,3) = ActionChargesList{nn}(:,3)-sum(repmat(ActionChargesList{nn}(:,3),[1,numel(DropChargesList)])>repmat(DropChargesList,[size(ActionChargesList{nn},1),1]),2);
                end
                
                ChargeExtTempC = ChargeExtTempC(MultList,:);
                
                for kk = 3:(size(ChargeExtTempC,2)-1)
                    [Temp, MultList, Multiplicities] = SymmetryUsed.FuseChargeList(ChargeIntTempC(:,kk-2), ChargeExtTempCAdjust(:,kk));
                    ChargeIntTempC = [ChargeIntTempC(MultList,:), Temp];
                    MultiplicitiesIntTempC = [MultiplicitiesIntTempC(MultList,:), Multiplicities];
                    
                    NumberList = 1:size(ChargeExtTempC,1);
                    DropChargesList = NumberList(~any(repmat(NumberList,[size(MultList,1),1])==repmat(MultList,[1,size(NumberList,2)]),1));
                    
                    if ~isempty(DropChargesList)
                        ActionChargesList{nn}(any(repmat(ActionChargesList{nn}(:,3),[1,numel(DropChargesList)])==repmat(DropChargesList,[size(ActionChargesList{nn},1),1]),2),:) = [];
                        %ActionChargesList{nn}(:,3) = ActionChargesList{nn}(:,3)-sum(repmat(ActionChargesList{nn}(:,3),[1,numel(DropChargesList)])>repmat(DropChargesList,[size(ActionChargesList{nn},1),1]),2);
                    end
                    
                    ChargeExtTempC = ChargeExtTempC(MultList,:);
                end
                %MultiplicitiesIntTempC = ones(size(ChargeIntTempC)+[0,1]);
            end
            
            Temp = ActionTotalTensors==ActionTotalTensors(nn,3);
            ActionTensorsChargesExternal(Temp) = repmat({ChargeExtTempC}, [sum(sum(Temp)),1]);
            ActionTensorsChargesInternal(Temp) = repmat({ChargeIntTempC}, [sum(sum(Temp)),1]);
            ActionTensorsMultiplicitiesInternal(Temp) = repmat({MultiplicitiesIntTempC}, [sum(sum(Temp)),1]);
        end
    else
        ActionTensorsConvertNumbers = cell(size(ActionTensorsConvertWords));
        for nn = 1:size(ActionTotalTensors,1)
            if ActionTypeTotal(nn) == 1 %trace
                LabelsA = [ActionTensorsChargesExternal{nn,1},ActionTensorsChargesInternal{nn,1},ActionTensorsMultiplicitiesInternal{nn,1}];
                
                ActionTensorsConvertNumbers{nn,3} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,3},...
                    ActionTensorsChargesExternal{nn,1},ActionTensorsChargesInternal{nn,1},ActionTensorsMultiplicitiesInternal{nn,1});
                
                %compute LabelsOut for trace:
                InLabels = ActionTensorsConvertNumbers{nn,3}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,3}{2};
                Moves = ActionTensorsConvertNumbers{nn,3}{3};
                
                
                TempLabels =[InLabels;LabelsA];
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows');
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOut = OutLabels(NumbersOut,:);
                
                [ListIn,ListOut,ListMult] = find(Moves(NumbersIn,NumbersOut));
                
                ActionChargesList{nn} = [ListIn, zeros(size(ListIn)),ListOut];
                ActionTensorsListMultNumbers{nn} = ListMult;
                
            elseif ActionTypeTotal(nn) == 2 %contraction
                LabelsA = [ActionTensorsChargesExternal{nn,1},ActionTensorsChargesInternal{nn,1},ActionTensorsMultiplicitiesInternal{nn,1}];
                LabelsB = [ActionTensorsChargesExternal{nn,2},ActionTensorsChargesInternal{nn,2},ActionTensorsMultiplicitiesInternal{nn,2}];
                
                %disp(['nn = ',num2str(nn)])%DebugDisp
                ActionTensorsConvertNumbers{nn,1} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,1},...
                    ActionTensorsChargesExternal{nn,1},ActionTensorsChargesInternal{nn,1},ActionTensorsMultiplicitiesInternal{nn,1});
                ActionTensorsConvertNumbers{nn,2} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,2},...
                    ActionTensorsChargesExternal{nn,2},ActionTensorsChargesInternal{nn,2},ActionTensorsMultiplicitiesInternal{nn,2});
                
                %Compute Reordered charges
                InLabels = ActionTensorsConvertNumbers{nn,1}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,1}{2};
                Moves = ActionTensorsConvertNumbers{nn,1}{3};
                
                
                TempLabels =[InLabels;LabelsA];
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOutA = OutLabels(NumbersOut,:);
                
                [ListInA,ListOutA,ListMultA] = find(Moves(NumbersIn,NumbersOut));
                
                %second one
                InLabels = ActionTensorsConvertNumbers{nn,2}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,2}{2};
                Moves = ActionTensorsConvertNumbers{nn,2}{3};
                
                
                TempLabels =[InLabels;LabelsB];
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOutB = OutLabels(NumbersOut,:);
                
                [ListInB,ListOutB,ListMultB] = find(Moves(NumbersIn,NumbersOut));
                
                
                %combine charges
                
                
                LegsNum = numel(ActualLegsTotal{nn,1});
                IntermediateAChargesExternal = LabelsOutA(:,1:LegsNum);
                if LegsNum >2
                    IntermediateAChargesInternal = LabelsOutA(:,LegsNum+(1:(LegsNum-2)));
                    IntermediateAMultiplicitiesInternal = LabelsOutA(:,2*LegsNum-2+(1:(LegsNum-2)));
                else
                    IntermediateAChargesInternal = zeros([size(LabelsOutA,1),0]);
                    IntermediateAMultiplicitiesInternal = zeros([size(LabelsOutA,1),0]);
                end
                
                
                LegsNum = numel(ActualLegsTotal{nn,2});
                IntermediateBChargesExternal = LabelsOutB(:,1:LegsNum);
                if LegsNum >2
                    IntermediateBChargesInternal = LabelsOutB(:,LegsNum+(1:(LegsNum-2)));
                    IntermediateBMultiplicitiesInternal = LabelsOutB(:,2*LegsNum-2+(1:(LegsNum-2)));
                else
                    IntermediateBChargesInternal = zeros([size(LabelsOutB,1),0]);
                    IntermediateBMultiplicitiesInternal = zeros([size(LabelsOutB,1),0]);
                end
                
                NumberLabelsAllA = size(IntermediateAChargesExternal,1);
                NumberLabelsAllB = size(IntermediateBChargesExternal,1);
                NumberExtLegsAllA = size(IntermediateAChargesExternal,2);
                NumberExtLegsAllB = size(IntermediateBChargesExternal,2);
                NumberIntLegsAllA = size(IntermediateAChargesInternal,2);
                NumberIntLegsAllB = size(IntermediateBChargesInternal,2);
                
                NumberMultsAllA = size(IntermediateAMultiplicitiesInternal,2);
                NumberMultsAllB = size(IntermediateBMultiplicitiesInternal,2);
                
                IntermediateAList = reshape(repmat(reshape(1:NumberLabelsAllA,[1,NumberLabelsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,1]);
                IntermediateAChargesExternal = reshape(repmat(reshape(IntermediateAChargesExternal,[1,NumberLabelsAllA,NumberExtLegsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,NumberExtLegsAllA]);
                IntermediateAChargesInternal = reshape(repmat(reshape(IntermediateAChargesInternal,[1,NumberLabelsAllA,NumberIntLegsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,NumberIntLegsAllA]);
                IntermediateAMultiplicitiesInternal = reshape(repmat(reshape(IntermediateAMultiplicitiesInternal,[1,NumberLabelsAllA,NumberMultsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,NumberMultsAllA]);
                
                IntermediateBList = reshape(repmat(reshape(1:NumberLabelsAllB,[NumberLabelsAllB,1]),...
                    [1,NumberLabelsAllA]),[NumberLabelsAllB*NumberLabelsAllA,1]);
                IntermediateBChargesExternal = reshape(repmat(reshape(IntermediateBChargesExternal,[NumberLabelsAllB,1,NumberExtLegsAllB]),...
                    [1,NumberLabelsAllA]),[NumberLabelsAllB*NumberLabelsAllA,NumberExtLegsAllB]);
                IntermediateBChargesInternal = reshape(repmat(reshape(IntermediateBChargesInternal,[NumberLabelsAllB,1,NumberIntLegsAllB]),...
                    [1,NumberLabelsAllA]),[NumberLabelsAllB*NumberLabelsAllA,NumberIntLegsAllB]);
                IntermediateBMultiplicitiesInternal = reshape(repmat(reshape(IntermediateBMultiplicitiesInternal,[NumberLabelsAllB,1,NumberMultsAllB]),...
                    [1,NumberLabelsAllA]),[NumberLabelsAllB*NumberLabelsAllA,NumberMultsAllB]);
                
                IntermediateABChargesExternal = [IntermediateAChargesExternal,IntermediateBChargesExternal];
                IntermediateABChargesInternal = [IntermediateAChargesInternal,IntermediateBChargesInternal];
                IntermediateABMultiplicitiesInternal = [IntermediateAMultiplicitiesInternal,IntermediateBMultiplicitiesInternal];
                
                clear IntermediateAChargesExternal IntermediateAChargesInternal IntermediateAMultiplicitiesInternal
                clear IntermediateBChargesExternal IntermediateBChargesInternal IntermediateBMultiplicitiesInternal
                
                IntermediateCChargesInternal = repmat(SymmetryUsed.getTrivialIrrep, ...
                    [size(IntermediateABChargesInternal,1), size(ActionTensorsCombineCharges{nn,2},2)]);
                IntermediateCMultiplicitiesInternal = ones([size(IntermediateABMultiplicitiesInternal,1), size(ActionTensorsCombineCharges{nn,3},2)]);
                
                IntermediateCChargesExternal = IntermediateABChargesExternal(:,ActionTensorsCombineCharges{nn,1});
                IntermediateCChargesInternal(:,ActionTensorsCombineCharges{nn,2}>0) = IntermediateABChargesInternal(:,ActionTensorsCombineCharges{nn,2}(ActionTensorsCombineCharges{nn,2}>0));
                IntermediateCChargesInternal(:,ActionTensorsCombineCharges{nn,2}<0) = IntermediateABChargesExternal(:,-ActionTensorsCombineCharges{nn,2}(ActionTensorsCombineCharges{nn,2}<0));
                IntermediateCMultiplicitiesInternal(:,ActionTensorsCombineCharges{nn,3}>0) = IntermediateABMultiplicitiesInternal(:,ActionTensorsCombineCharges{nn,3}(ActionTensorsCombineCharges{nn,3}>0));
                
                
                clear IntermediateABChargesExternal IntermediateABChargesInternal IntermediateABMultiplicitiesInternal
                
                %this stores everything and therefore is already unique.
                
                %ActionTensorsConvertNumbers{nn,3} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,3},...
                %       IntermediateCChargesExternal,IntermediateCChargesInternal,IntermediateCMultiplicitiesInternal);
                
                KeptSame = ActionTensorsConvertWords{nn,3}{1,2};
                Keep = true([size(IntermediateCChargesExternal,1),1]);
                for aa = 1:size(KeptSame,2)
                    if KeptSame(3,aa)*KeptSame(4,aa)==+1
                        Keep = Keep&(IntermediateCChargesExternal(:,-KeptSame(1,aa))==IntermediateCChargesExternal(:,-KeptSame(2,aa)));
                    else
                        Keep = Keep&(IntermediateCChargesExternal(:,-KeptSame(1,aa))==InverseIrrep(IntermediateCChargesExternal(:,-KeptSame(2,aa))));
                    end
                end
                
                IntermediateCChargesExternal = IntermediateCChargesExternal(Keep,:);
                IntermediateCChargesInternal = IntermediateCChargesInternal(Keep,:);
                IntermediateCMultiplicitiesInternal = IntermediateCMultiplicitiesInternal(Keep,:);
                IntermediateAList = IntermediateAList(Keep,:);
                IntermediateBList = IntermediateBList(Keep,:);
                
                %%{
                if nn ~=1;
 %                   toc;
                end
                
%                tic
                if size(IntermediateCChargesExternal,1)<=MaxChoices
                    ActionTensorsConvertNumbers{nn,3} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,3}(2:end,:),...
                       IntermediateCChargesExternal,IntermediateCChargesInternal,IntermediateCMultiplicitiesInternal);
                else
                    
                    Temp = cell([1,ceil(size(IntermediateCChargesExternal,1)/MaxChoices)]);
                    WorkingFirstLabels = [];
                    WorkingSecondLabels = [];
                    WorkingMoves = [];
                    
                    for kk = 1:ceil(size(IntermediateCChargesExternal,1)/MaxChoices)
                        Working = (MaxChoices*(kk-1)+1):min(MaxChoices*kk,size(IntermediateCChargesExternal,1));
                        
                        Temp = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,3},...
                        IntermediateCChargesExternal(Working,:),IntermediateCChargesInternal(Working,:),IntermediateCMultiplicitiesInternal(Working,:));
                        if kk == 1
                            WorkingFirstLabels = Temp{1};
                            WorkingSecondLabels = Temp{2};
                            WorkingMoves = Temp{3};
                        else
                            SizesFirst = size(WorkingMoves,1);
                            SizesSecond = size(WorkingMoves,2);
                            
                            [WorkingFirstLabels, ~, AdjustFirst] = unique([WorkingFirstLabels;Temp{1}],'rows');
                            [WorkingSecondLabels, ~, AdjustSecond] = unique([WorkingSecondLabels;Temp{2}],'rows');
                            
                            WorkingAdjustFirst = AdjustFirst(1:SizesFirst);
                            WorkingAdjustSecond = AdjustSecond(1:SizesSecond);
                            
                            TempAdjustFirst = AdjustFirst((SizesFirst+1):end);
                            TempAdjustSecond = AdjustSecond((SizesSecond+1):end);
                            
                            TempMatrix = sparse(size(WorkingFirstLabels,1),size(WorkingSecondLabels,1));
                            TempWorkingMatrix = sparse(size(WorkingFirstLabels,1),size(WorkingSecondLabels,1));
                            
                            TempWorkingMatrix(WorkingAdjustFirst,WorkingAdjustSecond) = WorkingMoves;
                            TempMatrix(TempAdjustFirst,TempAdjustSecond) = Temp{3};
                            
                            WorkingMoves = TempMatrix+TempWorkingMatrix;
                            
                        end
                    end
                    
                    ActionTensorsConvertNumbers{nn,3} = {WorkingFirstLabels, WorkingSecondLabels, WorkingMoves};
                    
                end
                %}
                
                
                %compute new Numbers for after combination
                
                %work out how the numbers then change:
                
                InLabels = ActionTensorsConvertNumbers{nn,3}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,3}{2};
                Moves = ActionTensorsConvertNumbers{nn,3}{3};
                
                
                TempLabels =[InLabels;[IntermediateCChargesExternal,IntermediateCChargesInternal,IntermediateCMultiplicitiesInternal]];
                
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
                
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOut = OutLabels(NumbersOut,:);
                
                [ListInC,ListOutC,ListMultC] = find(Moves(NumbersIn,NumbersOut));
                
                %now combine the ListMults:
                ListInterA = IntermediateAList(ListInC);%ListOutA(IntermediateAList(ListInC));
                ListInterB = IntermediateBList(ListInC);%ListOutB(IntermediateBList(ListInC));
                
                %now need to factor in the initial multiplications that was done to A and B (using ListInA,ListOutA,ListMultA and ListInB,ListOutB,ListMultB).
                
                MaxSparseAIn = max(ListInA);
                MaxSparseAOut = max(ListOutA);
                
                MaxSparseBIn = max(ListInB);
                MaxSparseBOut = max(ListOutB);
                
                MaxSparseCOut = max(ListOutC);
                
                
                SparseA = sparse(ListInA, ListOutA, ListMultA,MaxSparseAIn,MaxSparseAOut);
                SparseB = sparse(ListInB, ListOutB, ListMultB,MaxSparseBIn,MaxSparseBOut);
                
                SparseC = sparse(ListInterA, (ListInterB-1)*MaxSparseCOut+ListOutC, ListMultC,MaxSparseAOut,MaxSparseBOut*MaxSparseCOut);
                [ListInModA,ListOutModBC,ListMultCMod] = find(SparseA*SparseC);
                ListOutModC = mod(ListOutModBC-1,MaxSparseCOut)+1;
                ListInModB = (ListOutModBC-ListOutModC)/MaxSparseCOut+1;
                if ~IsInteger(ListInModB)
                    error('Affirmation Error: ListInModB should be integers');
                end
                
                SparseCMod = sparse(ListInModB, (ListInModA-1)*MaxSparseCOut+ListOutModC, ListMultCMod,MaxSparseBOut,MaxSparseAIn*MaxSparseCOut);
                [ListInFinalB,ListOutFinalAC,ListMultFinal] = find(SparseB*SparseCMod);
                ListOutFinalC = mod(ListOutFinalAC-1,MaxSparseCOut)+1;
                ListInFinalA = (ListOutFinalAC-ListOutFinalC)/MaxSparseCOut+1;
                if ~IsInteger(ListInFinalA)
                    error('Affirmation Error: ListInModB should be integers');
                end
                %{
AAPrime = full(ActionTensorsConvertNumbers{nn,1}{3});
BBPrime = full(ActionTensorsConvertNumbers{nn,2}{3});
AANum1 = size(ActionTensorsConvertNumbers{nn,1}{1},1);
BBNum1 = size(ActionTensorsConvertNumbers{nn,2}{1},1);
AANum2 = size(ActionTensorsConvertNumbers{nn,1}{2},1);
BBNum2 = size(ActionTensorsConvertNumbers{nn,2}{2},1);
ABABNum = size(ActionTensorsConvertNumbers{nn,3}{1},1);
CCNum = size(ActionTensorsConvertNumbers{nn,3}{2},1);
if ABABNum~=AANum2*BBNum2; disp('error: A and B numbers don''t line up'); end
CCPrime = permute(reshape(full(ActionTensorsConvertNumbers{nn,3}{3}),[BBNum2,AANum2,CCNum]),[2,1,3]);
CC2Prime = reshape(AAPrime*reshape(CCPrime,[AANum2,BBNum2*CCNum]),[AANum1,BBNum2,CCNum]);
CC3Prime = permute(reshape(BBPrime*reshape(permute(CC2Prime,[2,1,3]),[BBNum2,AANum1*CCNum]),[BBNum1,AANum1,CCNum]),[2,1,3]);
SparseCFinal = sparse(ListInFinalA, (ListInFinalB-1)*MaxSparseCOut+ListOutFinalC, ListMultFinal,MaxSparseAIn,MaxSparseBIn*MaxSparseCOut);
CCFinal = permute(reshape(full(SparseCFinal),[AANum1,CCNum,BBNum1]),[1,3,2]);

[ListInFinalAB,ListOutFinalC,ListMultFinal] = find(reshape(CCFinal,[AANum1*BBNum1,CCNum]));

ListInFinalA = mod(ListInFinalAB-1,AANum1)+1;
ListInFinalB = (ListInFinalAB-ListInFinalA)/AANum1+1;
                if ~IsInteger(ListInFinalB)
                    error('Affirmation Error: ListInFinalB should be integers');
                end
                %}
                
                ActionChargesList{nn} = [ListInFinalA(:), ListInFinalB(:),ListOutFinalC(:)];
                ActionTensorsListMultNumbers{nn} = ListMultFinal;
  %              toc
                
            else %kronica
                LabelsA = [ActionTensorsChargesExternal{nn,1},ActionTensorsChargesInternal{nn,1},ActionTensorsMultiplicitiesInternal{nn,1}];
                LabelsB = [ActionTensorsChargesExternal{nn,2},ActionTensorsChargesInternal{nn,2},ActionTensorsMultiplicitiesInternal{nn,2}];
                
                ActionTensorsConvertNumbers{nn,1} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,1},...
                    ActionTensorsChargesExternal{nn,1},ActionTensorsChargesInternal{nn,1},ActionTensorsMultiplicitiesInternal{nn,1});
                ActionTensorsConvertNumbers{nn,2} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,2},...
                    ActionTensorsChargesExternal{nn,2},ActionTensorsChargesInternal{nn,2},ActionTensorsMultiplicitiesInternal{nn,2});
                
                %Compute Reordered charges
                InLabels = ActionTensorsConvertNumbers{nn,1}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,1}{2};
                Moves = ActionTensorsConvertNumbers{nn,1}{3};
                
                
                TempLabels =[InLabels;LabelsA];
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOutA = OutLabels(NumbersOut,:);
                
                [ListInA,ListOutA,ListMultA] = find(Moves(NumbersIn,NumbersOut));
                
                %second one
                InLabels = ActionTensorsConvertNumbers{nn,2}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,2}{2};
                Moves = ActionTensorsConvertNumbers{nn,2}{3};
                
                
                TempLabels =[InLabels;LabelsB];
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOutB = OutLabels(NumbersOut,:);
                
                [ListInB,ListOutB,ListMultB] = find(Moves(NumbersIn,NumbersOut));
                
                
                %combine charges
                
                
                LegsNum = numel(ActualLegsTotal{nn,1});
                IntermediateAChargesExternal = LabelsOutA(:,1:LegsNum);
                if LegsNum >2
                    IntermediateAChargesInternal = LabelsOutA(:,LegsNum+(1:(LegsNum-2)));
                    IntermediateAMultiplicitiesInternal = LabelsOutA(:,2*LegsNum-2+(1:(LegsNum-2)));
                else
                    IntermediateAChargesInternal = zeros([size(LabelsOutA,1),0]);
                    IntermediateAMultiplicitiesInternal = zeros([size(LabelsOutA,1),0]);
                end
                
                
                LegsNum = numel(ActualLegsTotal{nn,2});
                IntermediateBChargesExternal = LabelsOutB(:,1:LegsNum);
                if LegsNum >2
                    IntermediateBChargesInternal = LabelsOutB(:,LegsNum+(1:(LegsNum-2)));
                    IntermediateBMultiplicitiesInternal = LabelsOutB(:,2*LegsNum-2+(1:(LegsNum-2)));
                else
                    IntermediateBChargesInternal = zeros([size(LabelsOutB,1),0]);
                    IntermediateBMultiplicitiesInternal = zeros([size(LabelsOutB,1),0]);
                end
                
                NumberLabelsAllA = size(IntermediateAChargesExternal,1);
                NumberLabelsAllB = size(IntermediateBChargesExternal,1);
                NumberExtLegsAllA = size(IntermediateAChargesExternal,2);
                NumberExtLegsAllB = size(IntermediateBChargesExternal,2);
                NumberIntLegsAllA = size(IntermediateAChargesInternal,2);
                NumberIntLegsAllB = size(IntermediateBChargesInternal,2);
                
                NumberMultsAllA = size(IntermediateAMultiplicitiesInternal,2);
                NumberMultsAllB = size(IntermediateBMultiplicitiesInternal,2);
                
                
                IntermediateAList = reshape(repmat(reshape(1:NumberLabelsAllA,[1,NumberLabelsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,1]);
                IntermediateAChargesExternal = reshape(repmat(reshape(IntermediateAChargesExternal,[1,NumberLabelsAllA,NumberExtLegsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,NumberExtLegsAllA]);
                IntermediateAChargesInternal = reshape(repmat(reshape(IntermediateAChargesInternal,[1,NumberLabelsAllA,NumberIntLegsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,NumberIntLegsAllA]);
                IntermediateAMultiplicitiesInternal = reshape(repmat(reshape(IntermediateAMultiplicitiesInternal,[1,NumberLabelsAllA,NumberMultsAllA]),...
                    [NumberLabelsAllB,1]),[NumberLabelsAllA*NumberLabelsAllB,NumberMultsAllA]);
                
                IntermediateBList = reshape(repmat(reshape(1:NumberLabelsAllB,[1,NumberLabelsAllB]),...
                    [NumberLabelsAllA,1]),[NumberLabelsAllB*NumberLabelsAllA,1]);
                IntermediateBChargesExternal = reshape(repmat(reshape(IntermediateBChargesExternal,[1,NumberLabelsAllB,NumberExtLegsAllB]),...
                    [NumberLabelsAllA,1]),[NumberLabelsAllB*NumberLabelsAllA,NumberExtLegsAllB]);
                IntermediateBChargesInternal = reshape(repmat(reshape(IntermediateBChargesInternal,[1,NumberLabelsAllB,NumberIntLegsAllB]),...
                    [NumberLabelsAllA,1]),[NumberLabelsAllB*NumberLabelsAllA,NumberIntLegsAllB]);
                IntermediateBMultiplicitiesInternal = reshape(repmat(reshape(IntermediateBMultiplicitiesInternal,[1,NumberLabelsAllB,NumberMultsAllB]),...
                    [NumberLabelsAllA,1]),[NumberLabelsAllB*NumberLabelsAllA,NumberMultsAllB]);
                
                IntermediateABChargesExternal = [IntermediateAChargesExternal,IntermediateBChargesExternal];
                IntermediateABChargesInternal = [IntermediateAChargesInternal,IntermediateBChargesInternal];
                IntermediateABMultiplicitiesInternal = [IntermediateAMultiplicitiesInternal,IntermediateBMultiplicitiesInternal];
                
                IntermediateCChargesInternal = repmat(SymmetryUsed.getTrivialIrrep, ...
                    [size(IntermediateABChargesInternal,1), size(ActionTensorsCombineCharges{nn,2},2)]);
                IntermediateCMultiplicitiesInternal = ones([size(IntermediateABMultiplicitiesInternal,1), size(ActionTensorsCombineCharges{nn,3},2)]);
                
                IntermediateCChargesExternal = IntermediateABChargesExternal(:,ActionTensorsCombineCharges{nn,1});
                IntermediateCChargesInternal(:,ActionTensorsCombineCharges{nn,2}>0) = IntermediateABChargesInternal(:,ActionTensorsCombineCharges{nn,2}(ActionTensorsCombineCharges{nn,2}>0));
                IntermediateCChargesInternal(:,ActionTensorsCombineCharges{nn,2}<0) = IntermediateABChargesExternal(:,-ActionTensorsCombineCharges{nn,2}(ActionTensorsCombineCharges{nn,2}<0));
                IntermediateCMultiplicitiesInternal(:,ActionTensorsCombineCharges{nn,3}>0) = IntermediateABMultiplicitiesInternal(:,ActionTensorsCombineCharges{nn,3}(ActionTensorsCombineCharges{nn,3}>0));
                
                %this stores everything and therefore is already unique.
                
                
                %compute new Numbers for after combination
                ActionTensorsConvertNumbers{nn,3} = SymmetryUsed.GenerateNumbers(ActionTensorsConvertWords{nn,3},...
                    IntermediateChargesExternal,IntermediateChargesInternal,IntermediateMultiplicitiesInternal);
                
                %work out how the numbers then change:
                
                InLabels = ActionTensorsConvertNumbers{nn,3}{1};
                OutLabels = ActionTensorsConvertNumbers{nn,3}{2};
                Moves = ActionTensorsConvertNumbers{nn,3}{3};
                
                
                TempLabels =[InLabels;[IntermediateCChargesExternal,IntermediateCChargesInternal,IntermediateCMultiplicitiesInternal]];
                [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
                if any(NumbersIn>size(InLabels,1))
                    error('Error: There are some unknown labels in here');
                end
                NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
            
                [~,NumbersOut,~] = find(Moves(NumbersIn,:));
                NumbersOut = unique(NumbersOut);
                LabelsOut = OutLabels(NumbersOut,:);
                
                [ListInC,ListOutC,ListMultC] = find(Moves(NumbersIn,NumbersOut));
                
                %now combine the ListMults:
                ListInterA = IntermediateAList(ListInC);%ListOutA(IntermediateAList(ListInC));
                ListInterB = IntermediateBList(ListInC);%ListOutB(IntermediateBList(ListInC));
                
                %now need to factor in the initial multiplications that was done to A and B (using ListInA,ListOutA,ListMultA and ListInB,ListOutB,ListMultB).
                
                MaxSparseAIn = max(ListInA);
                MaxSparseAOut = max(ListOutA);
                
                MaxSparseBIn = max(ListInB);
                MaxSparseBOut = max(ListOutB);
                
                MaxSparseCOut = max(ListOutC);
                
                
                SparseA = sparse(ListInA, ListOutA, ListMultA,MaxSparseAIn,MaxSparseAOut);
                SparseB = sparse(ListInB, ListOutB, ListMultB,MaxSparseBIn,MaxSparseBOut);
                
                SparseC = sparse(ListInterA, (ListInterB-1)*MaxSparseCOut+ListOutC, ListMultC,MaxSparseAOut,MaxSparseBOut*MaxSparseCOut);
                [ListInModA,ListOutModBC,ListMultCMod] = find(SparseA*SparseC);
                ListOutModC = mod(ListOutModBC-1,MaxSparseCOut)+1;
                ListInModB = (ListOutModBC-ListOutModC)/MaxSparseCOut+1;
                if ~IsInteger(ListInModB)
                    error('Affirmation Error: ListInModB should be integers');
                end
                
                SparseCMod = sparse(ListInModB, (ListInModA-1)*MaxSparseCOut+ListOutModC, ListMultCMod,MaxSparseBOut,MaxSparseAIn*MaxSparseCOut);
                [ListInFinalB,ListOutFinalAC,ListMultFinal] = find(SparseB*SparseCMod);
                ListOutFinalC = mod(ListOutFinalAC-1,MaxSparseCOut)+1;
                ListInFinalA = (ListOutFinalAC-ListOutFinalC)/MaxSparseCOut+1;
                if ~IsInteger(ListInFinalA)
                    error('Affirmation Error: ListInModB should be integers');
                end
                
                
                ActionChargesList{nn} = [ListInFinalA, ListInFinalB,ListOutFinalC];
                ActionTensorsListMultNumbers{nn} = ListMultFinal;
                
            end
            
            LegsNum = numel(ActualLegsTotal{nn,3});
            ChargeExtTempC = LabelsOut(:,1:LegsNum);
            if LegsNum >2
                if size(LabelsOut,2) ~= LegsNum*3-4
                    error('Affirmation Error: LegsNum does not agree with the LabelsOut length')
                end
                ChargeIntTempC = LabelsOut(:,LegsNum+(1:(LegsNum-2)));
                MultiplicitiesIntTempC = LabelsOut(:,2*LegsNum-2+(1:(LegsNum-2)));
            else
                if size(LabelsOut,2) ~= LegsNum
                    error('Affirmation Error: LegsNum does not agree with the LabelsOut length')
                end
                ChargeIntTempC = zeros([size(LabelsOut,1),0]);
                MultiplicitiesIntTempC = zeros([size(LabelsOut,1),0]);
            end
            
            Temp = ActionTotalTensors==ActionTotalTensors(nn,3);
            ActionTensorsChargesExternal(Temp) = repmat({ChargeExtTempC}, [sum(sum(Temp)),1]);
            ActionTensorsChargesInternal(Temp) = repmat({ChargeIntTempC}, [sum(sum(Temp)),1]);
            ActionTensorsMultiplicitiesInternal(Temp) = repmat({MultiplicitiesIntTempC}, [sum(sum(Temp)),1]);
        end
    end
    for nn = 1:size(ActionTotalTensors,1)
        [~,Index,~] = unique(ActionChargesList{nn}(:,3));
        C_List{nn} = ActionChargesList{nn}(Index,:);
        ActionMaxNumTensorParts(nn) = max(C_List{nn}(:,3));
    end
            
    
    if size(LegOrder,2)~=3
        LegOrder(:,3) = cell([size(LegOrder,1),1]);
    end
    
            
            
            
    OldTrueChargesExternal = cell([size(OldTrue,1),1]);
    OldTrueChargesInternal = cell([size(OldTrue,1),1]);
    OldTrueMultiplicitiesInternal = cell([size(OldTrue,1),1]);
    
    %I am assuming here that I can use Structure to work out the
    %contraction order while going left to right.
    
    for nn = 1:size(OldTrue,2)
        Temp = ActionTotalTensors(:,3) == OldTrue(nn);
        if sum(sum(Temp))~=1
            error('Affirmation Failed: There should only be one copy of this produced')
        end
        
        if ~SymmetryUsed.IsBraiding && ~SymmetryUsed.IsNonAbelian && FlagSpeedUpAbelian
            ChargeExtTempC = ActionTensorsChargesExternal{Temp,3}(:,OldTruePermute{nn});
            
            %adjust for changing
            InverseIrrepUsed = SymmetryUsed.InverseIrrep;
            
            if ~isempty(ChargeExtTempC)
                [~,Index] = sort(OldTrueStructure{nn}(OldTrueStructure{nn}<0),'descend');
                ChargeExtDirTempC = OldTrueChargeDirections{nn}(OldTrueStructure{nn}<0);
                ChargeExtDirTempC = reshape(ChargeExtDirTempC(Index),[1,numel(Index)]);
                
                SourceStructure = ActionTensorsStructure{Temp,3};
                [~,Index] = sort(SourceStructure(SourceStructure<0),'descend');
                
                ChargeExtDirTempOrigC = C_ChargeDirections{Temp}(SourceStructure<0);
                ChargeExtDirTempOrigC = reshape(ChargeExtDirTempOrigC(Index),[1,numel(Index)]);
                ChargeExtDirTempOrigC = reshape(ChargeExtDirTempOrigC(OldTruePermute{nn}),[1,numel(Index)]);
                
                %if they aren't the same going in as coming out:
                ChargeExtTempC(repmat(ChargeExtDirTempOrigC.*ChargeExtDirTempC == -1,[size(ChargeExtTempC,1),1])) = ...
                    InverseIrrepUsed(ChargeExtTempC(repmat(ChargeExtDirTempOrigC.*ChargeExtDirTempC == -1,[size(ChargeExtTempC,1),1])));
                
                
                [~,Index] = sort(OldTrueStructure{nn}(OldTrueStructure{nn}>0),'ascend');
                ChargeIntDirTempC = OldTrueChargeDirections{nn}(OldTrueStructure{nn}>0);
                ChargeIntDirTempC = reshape(ChargeIntDirTempC(Index),[1,numel(Index)]);
                
                
            end
            
            if size(ChargeExtTempC,2) == 0
                %then there is only one output (the numeric output)
                ChargeIntTempC = zeros([1,0]);
                MultiplictiesIntTempC = zeros([1,0]);
                if size(ChargeExtTempC,1) ~=1
                    error('Affirmation Failed: ChargeExtTempC is not one entry')
                end
                
            elseif size(ChargeExtTempC,2) == 1
                ChargeIntTempC = zeros([size(ChargeExtTempC,1),0]);
                MultiplictiesIntTempC = zeros([size(ChargeExtTempC,1),0]);
                %check if they are on the same side
                
                if prod(OldTrueChargeSides{nn}) == 1 %same side
                    if (any(ChargeExtTempC(:,1) ~= InverseIrrep(ChargeExtTempC(:,2)))&&prod(ChargeExtDirTempC)==1) || ...
                            (any(ChargeExtTempC(:,1) ~= ChargeExtTempC(:,2))&&prod(ChargeExtDirTempC)==-1)
                        error('Affirmation Failed: ChargeExtTempC doesn''t have Inverse Irreps')
                    end
                else
                    if (any(ChargeExtTempC(:,1) ~= InverseIrrep(ChargeExtTempC(:,2)))&&prod(ChargeExtDirTempC)==-1) || ...
                            (any(ChargeExtTempC(:,1) ~= ChargeExtTempC(:,2))&&prod(ChargeExtDirTempC)==+1)
                        error('Affirmation Failed: ChargeExtTempC doesn''t have Inverse Irreps')
                    end
                end
                
            else
                ChargeIntTempC = zeros([size(ChargeExtTempC,1),size(ChargeExtTempC,2)-2]);
                MultiplicitiesIntTempC = ones([size(ChargeExtTempC,1),size(ChargeExtTempC,2)-2]);
                StructureUsed = OldTrueStructure{nn};
                
                
                %work out Order
                OrderInternalLegs = [];
                StructureUsedTemp  = StructureUsed;
                
                while numel(OrderInternalLegs)<size(StructureUsed,2)
                    
                    Using = find(all(StructureUsedTemp<0,1));
                    
                    if isempty(Using)
                        error('Affirmation Error: Using is empty')
                    end
                    
                    OrderInternalLegs = [OrderInternalLegs,Using];
                    StructureUsedTemp(:,OrderInternalLegs) = 0;
                    for aa = Using
                        if aa>0
                            StructureUsedTemp(StructureUsedTemp == aa-1) = -aa;
                        end
                    end
                    
                end
                
                OrderInternalLegs = OrderInternalLegs-1;
                OrderInternalLegs(OrderInternalLegs == 0) = [];
                
                for kk = OrderInternalLegs+1%2:size(StructureUsed,2)
                    
                    if StructureUsed(1,kk)<0
                        ChargeIn1 = ChargeExtTempC(:,-StructureUsed(1,kk));
                        ChargeSide1 = OldTrueChargeSides{nn}(-StructureUsed(1,kk));
                        ChargeDirection1 = ChargeExtDirTempC(-StructureUsed(1,kk));
                    else
                        ChargeIn1 = ChargeIntTempC(:,StructureUsed(1,kk));
                        ChargeSide1 = OldTrueChargeSidesInt{nn}(StructureUsed(1,kk));
                        ChargeDirection1 = ChargeIntDirTempC(StructureUsed(1,kk));
                    end
                    if StructureUsed(2,kk)<0
                        ChargeIn2 = ChargeExtTempC(:,-StructureUsed(2,kk));
                        ChargeSide2 = OldTrueChargeSides{nn}(-StructureUsed(2,kk));
                        ChargeDirection2 = ChargeExtDirTempC(-StructureUsed(2,kk));
                    else
                        ChargeIn2 = ChargeIntTempC(:,StructureUsed(2,kk));
                        ChargeSide2 = OldTrueChargeSidesInt{nn}(StructureUsed(2,kk));
                        ChargeDirection2 = ChargeIntDirTempC(StructureUsed(2,kk));
                    end
                    
                    if kk ~= 1
                        ChargeSideOut = -OldTrueChargeSidesInt{nn}(kk-1);
                        ChargeDirectionOut = ChargeIntDirTempC(kk-1);
                        
                        
                        
                        
                        InverseDir1 = ChargeDirection1*ChargeDirectionOut;%going into the center is -1
                        InverseDir2 = ChargeDirection2*ChargeDirectionOut;%going into the center is -1
                        
                        
                        %now we assume that +1 means same direction, now
                        %adjust for the sides, if 1 and 2 are on the same
                        %side then this is fine (assuming out is the
                        %opposite, otherwise one of them needs to be
                        %inverted again
                        
                        if ChargeSide1 == ChargeSide2
                            
                            if ChargeSideOut  == ChargeSide1
                                error('Affirmation Error: OldTrue requires all the legs to be on the same side for some vertex.');
                            end
                            
                        else
                            
                            if ChargeSide1 == ChargeSideOut
                                InverseDir1 = -InverseDir1;
                            elseif ChargeSide2 == ChargeSideOut
                                InverseDir2 = -InverseDir2;
                            else
                                error('Affirmation Error: ChargeSide Out is wrong');
                            end
                            
                        end
                        
                        
                        InverseDir1 = InverseDir1 == -1;
                        InverseDir2 = InverseDir2 == -1;
                        
                        
                        
                        
                        if ~InverseDir1 && ~InverseDir2
                            ChargeIntTempC(:,kk-1) = SymmetryUsed.FuseChargeList(ChargeIn1,ChargeIn2);
                        elseif ~InverseDir1 
                            ChargeIntTempC(:,kk-1) = SymmetryUsed.FuseChargeList(ChargeIn1,reshape(InverseIrrepUsed(ChargeIn2),size(ChargeIn2)));
                        elseif ~InverseDir2
                            ChargeIntTempC(:,kk-1) = SymmetryUsed.FuseChargeList(reshape(InverseIrrepUsed(ChargeIn1),size(ChargeIn2)),ChargeIn2);
                        else
                            ChargeIntTempC(:,kk-1) = SymmetryUsed.FuseChargeList(reshape(InverseIrrepUsed(ChargeIn1),size(ChargeIn2)),...
                                reshape(InverseIrrepUsed(ChargeIn2),size(ChargeIn2)));
                        end
                        
                        
                    else %kk == 1
                        ChargeDirectionCombo = ChargeDirection1 * ChargeDirection2;
                        
                        if ChargeSide1*ChargeSide2 == 1 %same side
                            if (any(ChargeIn1 ~= InverseIrrep(ChargeIn2))&&ChargeDirectionCombo==1) || ...
                                    (any(ChargeIn1 ~= ChargeIn2)&&ChargeDirectionCombo==-1)
                                error('Affirmation Failed: ChargeExtTempC doesn''t have Inverse Irreps')
                            end
                        else
                            if (any(ChargeIn1 ~= InverseIrrep(ChargeIn2))&&ChargeDirectionCombo==-1) || ...
                                    (any(ChargeIn1 ~= ChargeIn2)&&ChargeDirectionCombo==+1)
                                error('Affirmation Failed: ChargeExtTempC doesn''t have Inverse Irreps')
                            end
                        end
                        
                        
                    end
                    
                    
                end
            end
        
        else
            ChargeExtTempA = ActionTensorsChargesExternal{Temp,3};%(:,OldTruePermute{nn});
            ChargeIntTempA = ActionTensorsChargesInternal{Temp,3};%(:,OldTruePermute{nn});
            MultIntTempA = ActionTensorsMultiplicitiesInternal{Temp,3};%(:,OldTruePermute{nn});
            
            OldTrueConvertNumbers{nn} = SymmetryUsed.GenerateNumbers(OldTrueConvertWords{nn},...
                ChargeExtTempA, ChargeIntTempA, MultIntTempA);
            
            
            LabelsA = [ChargeExtTempA, ChargeIntTempA, MultIntTempA];
            
            %compute LabelsOut for trace:
            InLabels = OldTrueConvertNumbers{nn}{1};
            OutLabels = OldTrueConvertNumbers{nn}{2};
            Moves = OldTrueConvertNumbers{nn}{3};
            
            
            TempLabels =[InLabels;LabelsA];
            [~,LocationIn,NumbersIn] = unique(TempLabels,'rows','first');
            if any(NumbersIn>size(InLabels,1))
                error('Error: There are some unknown labels in here');
            end
            NumbersIn = LocationIn(NumbersIn((size(InLabels,1)+1):end));
            [~,NumbersOut,~] = find(Moves(NumbersIn,:));
            NumbersOut = unique(NumbersOut);
            %[~,NumbersOut,~] = unique(NumbersOut);
            LabelsOut = OutLabels(NumbersOut,:);
            
            [ListIn,ListOut,ListMult] = find(Moves(NumbersIn,NumbersOut));
            
            OldTrueList{nn} = [ListIn, zeros(size(ListIn)),ListOut];
            OldTrueListMultNumbers{nn} = ListMult;
            LegsNum = numel(ActualLegsTotal{Temp,3});
            ChargeExtTempC = LabelsOut(:,1:LegsNum);
            if LegsNum >2
                if size(LabelsOut,2) ~= LegsNum*3-4
                    error('Affirmation Error: LegsNum does not agree with the LabelsOut length')
                end
                ChargeIntTempC = LabelsOut(:,LegsNum+(1:(LegsNum-2)));
                MultiplicitiesIntTempC = LabelsOut(:,2*LegsNum-2+(1:(LegsNum-2)));
            else
                if size(LabelsOut,2) ~= LegsNum
                    error('Affirmation Error: LegsNum does not agree with the LabelsOut length')
                end
                ChargeIntTempC = zeros([size(LabelsOut,1),0]);
                MultiplicitiesIntTempC = zeros([size(LabelsOut,1),0]);
            end
            
        end
        
        OldTrueChargesExternal{nn} = ChargeExtTempC;
        OldTrueChargesInternal{nn} = ChargeIntTempC;
        OldTrueMultiplicitiesInternal{nn} = MultiplicitiesIntTempC;
        
    end
    
    
    if ~isequal(SymmetryUsed.getSymName,'none')
        C_MultiplicitiesInternal = ActionTensorsMultiplicitiesInternal(:,3);
        C_ChargeLabelsInternal = ActionTensorsChargesInternal(:,3);
        C_ChargeLabelsExternal = ActionTensorsChargesExternal(:,3);
    end
    %SavedInfoMidDetails = [C_Structure, C_ChargeDirections, C_MultiplicitiesInternal, C_ChargeLabelsInternal, C_ChargeLabelsExternal, ActionChargesList, C_List, ActionMaxNumTensorParts,...
%                                  OldTrueStructure,OldTrueChargeDimensions,OldTrueChargesExternal,OldTrueChargesInternal,OldTrueMultiplicitiesInternal];
    
end

if FlagLoadSavePlanarMid
    
    %this is just a repeat of the tensor thing in PlanarMid
    [OldTrueConvertNumbers, ActionTensorsConvertWords, ActionTensorsCombineCharges, OldTrueListMultNumbers, ActionTensorsListMultNumbers,...
                OldTrueList,ActionChargesList,ActionTensorsChargeDirectionsExternalOriginal] = SymmetryUsed.LoadCurrentPlanarMidData(FlagNumber);
end

if FlagLoadSaveMid
    [C_MultiplicitiesInternal, C_ChargeLabelsInternal, C_ChargeLabelsExternal, ActionChargesList, C_List, ActionMaxNumTensorParts,...
                                  OldTrueChargesExternal,OldTrueChargesInternal,OldTrueMultiplicitiesInternal] = SymmetryUsed.LoadCurrentMidData(FlagNumber);
end


if FlagNewFineDetails
    
    
    if FlagLoadSaveBroad
        [UniqueLegsPos] = SymmetryUsed.LoadDataForFineWithoutBroad(FlagNumber);
    else
        if ~isempty(UniqueLegsNeg) && FlagNewPlanarDetails
            UniqueLegsPos = UniqueLegsPosOld;
        end
    end
    
    for kk = UniqueLegsPos
    
    %first find the two legs:
    
    NumberFound = 0;
    Tensor1 = 0;
    Tensor2 = 0;
    Leg1 = 0;
    Leg2 = 0;
    
    for ll = 1:length(Legs)
        TempValues = Legs{ll} == kk;
        TempSum = sum(sum(TempValues));
        if TempSum == 1;
            if NumberFound == 0
                Tensor1 = ll;
                Leg1 = find(any(TempValues,1));
                %Note that the copy doesn't matter due to the fact that all of those tensors are the same 
                NumberFound = 1;
                
            else %then NumberFound == 1
                Tensor2 = ll;
                Leg2 = find(any(TempValues,1));
                %Note that the copy doesn't matter due to the fact that all of those tensors are the same
                NumberFound = 2;
                
            end
            
            
            
            
        elseif TempSum == 2;
            
            %then we haven't found any copies so far.
            Tensor1 = ll;
            Tensor2 = ll;
            
            LegLocations2 = find(any(TempValues,1));
            
            if numel(LegLocations2) == 2;
                Leg1 = LegLocations2(1);
                Leg2 = LegLocations2(2);
            else %numel(LegLocations2) = 1
                Leg1 = LegLocations2(1);
                Leg2 = LegLocations2(1);
            end
            
            %Note that the copy doesn't matter due to the fact that all of those tensors are the same
            NumberFound = 2;
        end
        
        if NumberFound == 2;
            
            if ~ConnectableLegs(Tensors{Tensor1}, Tensors{Tensor2}, Leg1, Leg2)
                error('We need the legs of Tensors %i and %i to agree on bond %i', Tensor1, Tensor2, kk);
            end
            
            break;
        end
        
    end
    
    if NumberFound ~= 2
        error('SYSTEM ERROR: Something went wrong when trying to make sure the number of copies of each bond is exactly 2 in all cases')
    end
    
end

end
    
if FlagLoadSaveFine
    %currently nothing
end




if FlagSaveAtEnd
    SavedInfoDetails = {ActionTensorsLocations, ActionTypeTotal, ActionTensorsPermute, ActionNumberContractedLegs,...
                        OldTrueLocations, OldTruePermute, OldTrueChargeSides, OldTrueChargeSidesInt,...
                        OldTrueStructure,OldTrueChargeDirections,ActionTotalTensors,DropTensor,OldTrue,OldTrueBraiding,OldTrueBraidingDirection,...
                        UniqueLegsNeg, FlagUserDidNotChooseLegOrder,TotalTensors, ActualLegsTotal, LegOrder,ActionTensorsStructure,...
                        ActionTensorsChargeDirections, UniqueLegsPos};
    SymmetryUsed.SaveCurrentLoad(FlagNumber, SavedCheckDetails, SavedInfoDetails);
    
    if ~isequal(SymmetryUsed.getSymName,'none')
        SavedInfoMidDetails = {C_MultiplicitiesInternal, C_ChargeLabelsInternal, C_ChargeLabelsExternal, ActionChargesList, C_List, ActionMaxNumTensorParts,...
                             OldTrueChargesExternal,OldTrueChargesInternal,OldTrueMultiplicitiesInternal};
        SavedInfoFineDetails = {[]};
        SymmetryUsed.SaveCurrentLoadMid(FlagNumber, SavedCheckMidDetails, SavedInfoMidDetails);
        SymmetryUsed.SaveCurrentLoadFine(FlagNumber, SavedCheckFineDetails, SavedInfoFineDetails);
    end
end

if FlagSavePlanarAtEnd&&~FlagSaveMidPlanarAtEnd
    SavedInfoPlanarDetails = {OldTrueConvertWords, ActionTensorsConvertWords, ActionTensorsCombineCharges,[],[],[],ActionChargesList,...
        ActionTensorsChargeDirectionsExternalOriginal};
    SymmetryUsed.SaveCurrentLoadPlanar(FlagNumber, {}, SavedInfoPlanarDetails);
end

if FlagSaveMidPlanarAtEnd
    SavedInfoPlanarDetails = {OldTrueConvertWords, ActionTensorsConvertWords, ActionTensorsCombineCharges, OldTrueListMultNumbers, ActionTensorsListMultNumbers,...
                OldTrueList,ActionChargesList,ActionTensorsChargeDirectionsExternalOriginal};
    SymmetryUsed.SaveCurrentLoadPlanar(FlagNumber, {}, SavedInfoPlanarDetails);
end

if FlagSaveMidAtEnd
    SavedInfoMidDetails = {C_MultiplicitiesInternal, C_ChargeLabelsInternal, C_ChargeLabelsExternal, ActionChargesList, C_List, ActionMaxNumTensorParts,...
                        OldTrueChargesExternal,OldTrueChargesInternal,OldTrueMultiplicitiesInternal};
    SavedInfoFineDetails = {[]};
    SymmetryUsed.SaveCurrentLoadMid(FlagNumber, SavedCheckMidDetails, SavedInfoMidDetails);
    SymmetryUsed.SaveCurrentLoadFine(FlagNumber, SavedCheckFineDetails, SavedInfoFineDetails);
end


if FlagSaveFineAtEnd
    SavedInfoFineDetails = {[]};
    SymmetryUsed.SaveCurrentLoadFine(FlagNumber, SavedCheckFineDetails, SavedInfoFineDetails);
end

%now the useful things are:
    % - ActionTensorsLocations: An Actions-by-3 array of where to look and
    %store the output tensors
    % - ActualLegsTotal: An actions-by-3 cell of the leg lables for these
    %tensors
    % - ActionTypesTotal: An Actions array of what type of contractions to do
    % - LegsContractedTotal: an Actions cell of the legs to contract for this
    %action
    % - KeepUnique: a Actions-by-2 list of if we should dump this site
    % - DropTensor: the negative of KeepUnique
    % - OldTrueLocations: a List of locations for where to look for each
    %environment output at the end of each procedure
    % - PermuteFinal: A List of permutations required for each final tensor
    % - ActionTensorsConvert: How to convert the tensors
    % - ActionTensorsPermute: How to permute the tensors
    % - ActionNumberContractedLegs: How many legs are being contracted
    % - OldTrueConvert: How to convert before outputting
    % - OldTruePermute: How to permute before outputting, identical to
    % PermuteFinal
    






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now the computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is based on entries in: 

TensorsRaw = cell(size(Tensors));
for kk = 1:numel(Tensors)
    TensorsRaw{kk} = Tensors{kk}.TensorEntries;
    TensorsSizes{kk} = Tensors{kk}.TensorEntriesSizes;
end

clear Tensors;

ActionTensorsForgottenLabels = cell(size(ActionTotalTensors));
OldTrueForgottenLabels = cell(size(OldTrue));
%note that there are no forgotten labels for the initials (this will need
%to be fixed later)








for nn = 1:size(ActionTensorsLocations,1)
    
    if ActionTypeTotal(nn) == 1
        %then we know that we are doing a trace
        
        %first take in what I'm going to be working with and perform the
        %moves so that the Tensor entries are in the correct order
        ATensors = TensorsRaw{ActionTensorsLocations(nn,1)};
        ASizes = TensorsSizes{ActionTensorsLocations(nn,1)};
        
        if ~isempty(ActionTensorsPermute{nn,1})
            for kk = 1:numel(ATensors)
                ATensors{kk} = permute(ATensors{kk},ActionTensorsPermute{nn,1});
            end
            ASizes = ASizes(:,ActionTensorsPermute{nn,1});
        end
        
        %delete the stuff I can forget:
        TensorsRaw(ActionTensorsLocations(nn,[DropTensor(nn,:),false])) = [];
        TensorsSizes(ActionTensorsLocations(nn,[DropTensor(nn,:),false])) = [];
        
        %now we need to play with the blocks
        ASizesReshape = [prod(ASizes(:,1:ActionNumberContractedLegs(nn)),2),...
            prod(ASizes(:,(1:ActionNumberContractedLegs(nn))+ActionNumberContractedLegs(nn)),2),...
            prod(ASizes(:,(2*ActionNumberContractedLegs(nn)+1):end),2)];
        
        
        if isequal(SymmetryUsed.getSymName,'none')
            ASizesTrace = 1:(ASizesReshape(1)+1):ASizesReshape(1)^2;
            
            CSizes = [ASizes((2*ActionNumberContractedLegs(nn)+1):end)];
            if numel(CSizes)==0
                CSizesUsed = [1,1];
            elseif numel(CSizes) == 1
                CSizesUsed = [CSizes,1];
            else
                CSizesUsed = CSizes;
            end
            
            CTensors = reshape(ATensors{1},[ASizesReshape(1)^2,ASizesReshape(3)]);
            CTensors = {reshape(sum(CTensors(ASizesTrace,:),1),CSizesUsed)};
            
        else
            OldForgottenLabelsA = ActionTensorsForgottenLabels{nn,1};
            TranslateA = 1:(numel(ATensors)+numel(OldForgottenLabelsA));
            Temp = TranslateA; Temp(OldForgottenLabelsA) = [];
            TranslateA(Temp) = 1:numel(ATensors);
            
            CSizes = ASizes(C_List{nn}(:,1),(2*ActionNumberContractedLegs(nn)+1):end);
            if numel(CSizes)==0
                CSizesUsed = ones([size(CSizes,1),2]);
            elseif numel(CSizes) == 1
                CSizesUsed = [CSizes,ones([size(CSizes,1),1])];
            else
                CSizesUsed = CSizes;
            end
            List = ActionChargesList{nn};
            
            CTensors = cell([ActionMaxNumTensorParts(nn),1]);
            
            if ~FlagMustUsePlanar
                for kk = 1:size(List,1)
                    if ~any(List(kk,1)==OldForgottenLabelsA)
                        ASizesTrace = 1:(ASizesReshape(TranslateA(List(kk,1)),1)+1):ASizesReshape(TranslateA(List(kk,1)),1)^2;
                        Temp = reshape(ATensors{TranslateA(List(kk,1)),1},[ASizesReshape(TranslateA(List(kk,1)),1)^2,ASizesReshape(TranslateA(List(kk,1)),3)]);
                        if isempty(CTensors{List(kk,3)})
                            CTensors{List(kk,3)} = reshape(sum(Temp(ASizesTrace,:),1),CSizesUsed(List(kk,3),:));
                        else
                            CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(sum(Temp(ASizesTrace,:),1),CSizesUsed(List(kk,3),:));
                        end
                    end
                end
            else
                ListMultNumbers = ActionTensorsListMultNumbers{nn};
                for kk = 1:size(List,1)
                    if ~any(List(kk,1) == OldForgottenLabelsA)
                        ASizesTrace = 1:(ASizesReshape(TranslateA(List(kk,1)),1)+1):ASizesReshape(TranslateA(List(kk,1)),1)^2;
                        Temp = reshape(ATensors{TranslateA(List(kk,1)),1},[ASizesReshape(TranslateA(List(kk,1)),1)^2,ASizesReshape(TranslateA(List(kk,1)),3)]);
                        if isempty(CTensors{List(kk,3)})
                            CTensors{List(kk,3)} = reshape(sum(Temp(ASizesTrace,:),1),CSizesUsed(List(kk,3),:))...
                                *ListMultNumbers(kk);
                        else
                            CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(sum(Temp(ASizesTrace,:),1),CSizesUsed(List(kk,3),:))...
                                *ListMultNumbers(kk);
                        end
                    end
                end
            end
            
            NewForgottenLabels = zeros([1,0]);
            for ll = 1:numel(CTensors)
                if isempty(CTensors{ll})
                    NewForgottenLabels = [NewForgottenLabels,ll];
                elseif max(abs(CTensors{ll}(:)))<MinNonZero
                    NewForgottenLabels = [NewForgottenLabels,ll];
                end
            end
            
            CTensors(NewForgottenLabels) = [];
            CSizesUsed(NewForgottenLabels,:) = [];
            
            Temp = ActionTotalTensors == ActionTotalTensors(nn,3);
            ActionTensorsForgottenLabels(Temp) = repmat({NewForgottenLabels},[sum(Temp(:)),1]);
            
            TempOld = OldTrue == ActionTotalTensors(nn,3);
            if sum(TempOld(:))~=0
                OldTrueForgottenLabels(TempOld) = repmat({NewForgottenLabels},[sum(TempOld(:)),1]);
            end
            
        end
        
        %store the output:
        
        if ~isempty(ActionTensorsPermute{nn,3})
            for kk = 1:numel(CTensors)
                CTensors{kk} = permute(CTensors{kk}, ActionTensorsPermute{nn,3});
            end
            CSizesUsed = CSizesUsed(:,ActionTensorsPermute{nn,3});
        end
        
        TensorsRaw{ActionTensorsLocations(nn,3)} = CTensors;
        TensorsSizes{ActionTensorsLocations(nn,3)} = CSizesUsed;
        
    elseif ActionTypeTotal(nn) == 2
        %then we know we are taking a product
        
        
        %first take in what I'm going to be working with and perform the
        %moves so that the Tensor entries are in the correct order
        ATensors = TensorsRaw{ActionTensorsLocations(nn,1)};
        ASizes = TensorsSizes{ActionTensorsLocations(nn,1)};
        
        if ~isempty(ActionTensorsPermute{nn,1})
            for kk = 1:numel(ATensors)
                ATensors{kk} = permute(ATensors{kk},ActionTensorsPermute{nn,1});
            end
            ASizes = ASizes(:,ActionTensorsPermute{nn,1});
        end
        BTensors = TensorsRaw{ActionTensorsLocations(nn,2)};
        BSizes = TensorsSizes{ActionTensorsLocations(nn,2)};
        
        if ~isempty(ActionTensorsPermute{nn,2})
            for kk = 1:numel(BTensors)
                BTensors{kk} = permute(BTensors{kk},ActionTensorsPermute{nn,2});
            end
            BSizes = BSizes(:,ActionTensorsPermute{nn,2});
        end
        
        %delete the stuff I can forget:
        TensorsRaw(ActionTensorsLocations(nn,[DropTensor(nn,:),false])) = [];
        TensorsSizes(ActionTensorsLocations(nn,[DropTensor(nn,:),false])) = [];
        
        
        
        %now we need to play with the blocks
        ASizesReshape = [prod(ASizes(:,1:(end-ActionNumberContractedLegs(nn))),2),prod(ASizes(:,(end-ActionNumberContractedLegs(nn)+1):end),2)];
        BSizesReshape = [prod(BSizes(:,1:ActionNumberContractedLegs(nn)),2),prod(BSizes(:,(ActionNumberContractedLegs(nn)+1):end),2)];
        
        
        if isequal(SymmetryUsed.getSymName,'none')
            
            CSizes = [ASizes(1:(end-ActionNumberContractedLegs(nn))),BSizes((ActionNumberContractedLegs(nn)+1):end)];
            if numel(CSizes)==0
                CSizesUsed = [1,1];
            elseif numel(CSizes) == 1
                CSizesUsed = [CSizes,1];
            else
                CSizesUsed = CSizes;
            end
            
            CTensors = {reshape(reshape(ATensors{1},ASizesReshape)*reshape(BTensors{1},BSizesReshape),CSizesUsed)};
            
            
        else
            OldForgottenLabelsA = ActionTensorsForgottenLabels{nn,1};
            TranslateA = 1:(numel(ATensors)+numel(OldForgottenLabelsA));
            Temp = TranslateA; Temp(OldForgottenLabelsA) = [];
            TranslateA(Temp) = 1:numel(ATensors);
            
            OldForgottenLabelsB = ActionTensorsForgottenLabels{nn,2};
            TranslateB = 1:(numel(BTensors)+numel(OldForgottenLabelsB));
            Temp = TranslateB; Temp(OldForgottenLabelsB) = [];
            TranslateB(Temp) = 1:numel(BTensors);
            
            List = ActionChargesList{nn};
            
            KeepA = ~any(repmat(List(:,1),[1,numel(OldForgottenLabelsA)])==repmat(OldForgottenLabelsA(:)',[size(List,1),1]),2);
            KeepB = ~any(repmat(List(:,2),[1,numel(OldForgottenLabelsB)])==repmat(OldForgottenLabelsB(:)',[size(List,1),1]),2);
            KeepC = KeepA&KeepB;
            
            CListMod = List(KeepC,:);
            [~,UseListForC,~] = unique(CListMod(:,3));
            
            CSizes = zeros([size(C_List{nn},1),size(ASizes,2)+size(BSizes,2)-2*ActionNumberContractedLegs(nn)]);
            CSizes(CListMod(UseListForC,3),:) = [ASizes(TranslateA(CListMod(UseListForC,1)),1:(end-ActionNumberContractedLegs(nn))),...
                                 BSizes(TranslateB(CListMod(UseListForC,2)),(ActionNumberContractedLegs(nn)+1):end)];
                             
            if numel(CSizes)==0
                CSizesUsed = ones([size(CSizes,1),2]);
            elseif numel(CSizes) == 1
                CSizesUsed = [CSizes,ones([size(CSizes,1),1])];
            else
                CSizesUsed = CSizes;
            end
            
            if nn == 1
                disp('')
            end
            
            %%{
            for kk = unique(TranslateA(List(:,1)))
                ATensors{kk} = reshape(ATensors{kk},ASizesReshape(kk,:));
                
            end
            
            for kk = unique(TranslateB(List(:,2)))
                BTensors{kk} = reshape(BTensors{kk},BSizesReshape(kk,:));
            end
            %}
            
            CTensors = cell([ActionMaxNumTensorParts(nn),1]);
            if ~FlagMustUsePlanar
                for kk = 1:size(List,1)
                    if ~any(List(kk,1) == OldForgottenLabelsA) && ~any(List(kk,2) == OldForgottenLabelsB)
                        if isempty(CTensors{List(kk,3)})
                            CTensors{List(kk,3)} = ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))};
                            %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                            %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:));
                        else
                            CTensors{List(kk,3)} = CTensors{List(kk,3)} + ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))};
                            %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                            %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:));
                        end
                    end
                end
            else
                ListMultNumbers = ActionTensorsListMultNumbers{nn};
                for kk = 1:size(List,1)
                    %{
                    if ~any(List(kk,1) == OldForgottenLabelsA) && ~any(List(kk,2) == OldForgottenLabelsB)
                        ANumel = prod(ASizesReshape(TranslateA(List(kk,1)),:));
                        BNumel = prod(BSizesReshape(TranslateB(List(kk,2)),:));
                        CNumel = prod(CSizesUsed(List(kk,3),:));
                        if isempty(CTensors{List(kk,3)})
                            if CNumel<=ANumel && CNumel<=BNumel
                                CTensors{List(kk,3)} = (ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                                %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:))*ListMultNumbers(kk);
                            elseif ANumel <= BNumel
                                CTensors{List(kk,3)} = (ATensors{TranslateA(List(kk,1))}*ListMultNumbers(kk))*BTensors{TranslateB(List(kk,2))};
                                %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))*ListMultNumbers(kk)...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:));
                            else
                                CTensors{List(kk,3)} = ATensors{TranslateA(List(kk,1))}*(BTensors{TranslateB(List(kk,2))}*ListMultNumbers(kk));
                                %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:))*ListMultNumbers(kk),CSizesUsed(List(kk,3),:));
                            end
                        else
                            if CNumel<=ANumel && CNumel<=BNumel
                                CTensors{List(kk,3)} = CTensors{List(kk,3)} + (ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                                %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:))*ListMultNumbers(kk);
                            elseif ANumel <= BNumel
                                CTensors{List(kk,3)} = CTensors{List(kk,3)} + (ATensors{TranslateA(List(kk,1))}*ListMultNumbers(kk))*BTensors{TranslateB(List(kk,2))};
                                %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))*ListMultNumbers(kk)...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:));
                            else
                                CTensors{List(kk,3)} = CTensors{List(kk,3)} + ATensors{TranslateA(List(kk,1))}*(BTensors{TranslateB(List(kk,2))}*ListMultNumbers(kk));
                                %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:))*ListMultNumbers(kk),CSizesUsed(List(kk,3),:));
                            end
%                            CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
%                                *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:))*ListMultNumbers(kk);
                        end
                    end
                    %}
                    
                    if isempty(CTensors{List(kk,3)})
                        CTensors{List(kk,3)} = (ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                                
                        %{
                            if CNumel<=ANumel && CNumel<=BNumel
                                CTensors{List(kk,3)} = (ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                                %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:))*ListMultNumbers(kk);
                            elseif ANumel <= BNumel
                                CTensors{List(kk,3)} = (ATensors{TranslateA(List(kk,1))}*ListMultNumbers(kk))*BTensors{TranslateB(List(kk,2))};
                                %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))*ListMultNumbers(kk)...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:));
                            else
                                CTensors{List(kk,3)} = ATensors{TranslateA(List(kk,1))}*(BTensors{TranslateB(List(kk,2))}*ListMultNumbers(kk));
                                %CTensors{List(kk,3)} = reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:))*ListMultNumbers(kk),CSizesUsed(List(kk,3),:));
                            end
                        %}
                    else
                        CTensors{List(kk,3)} = CTensors{List(kk,3)} + (ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                                %{
                            if CNumel<=ANumel && CNumel<=BNumel
                                CTensors{List(kk,3)} = CTensors{List(kk,3)} + (ATensors{TranslateA(List(kk,1))}*BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                                %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:))*ListMultNumbers(kk);
                            elseif ANumel <= BNumel
                                CTensors{List(kk,3)} = CTensors{List(kk,3)} + (ATensors{TranslateA(List(kk,1))}*ListMultNumbers(kk))*BTensors{TranslateB(List(kk,2))};
                                %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))*ListMultNumbers(kk)...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:));
                            else
                                CTensors{List(kk,3)} = CTensors{List(kk,3)} + ATensors{TranslateA(List(kk,1))}*(BTensors{TranslateB(List(kk,2))}*ListMultNumbers(kk));
                                %CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
                                %    *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:))*ListMultNumbers(kk),CSizesUsed(List(kk,3),:));
                            end
%                            CTensors{List(kk,3)} = CTensors{List(kk,3)} + reshape(reshape(ATensors{TranslateA(List(kk,1))},ASizesReshape(TranslateA(List(kk,1)),:))...
%                                *reshape(BTensors{TranslateB(List(kk,2))},BSizesReshape(TranslateB(List(kk,2)),:)),CSizesUsed(List(kk,3),:))*ListMultNumbers(kk);
                        %}
                    end
                    
                end
            end
            
            %%{
            for kk = unique(List(:,3))'
                CTensors{kk} = reshape(CTensors{kk},CSizesUsed(kk,:));
            end
            %}
            
            NewForgottenLabels = zeros([1,0]);
            
            if SymmetryUsed.FlagSymCon_CheckIfEmpty
                for ll = 1:numel(CTensors)
                    if isempty(CTensors{ll})
                        NewForgottenLabels = [NewForgottenLabels,ll];
                    elseif max(abs(CTensors{ll}(:)))<MinNonZero
                        NewForgottenLabels = [NewForgottenLabels,ll];
                    end
                end
                CTensors(NewForgottenLabels) = [];
                CSizesUsed(NewForgottenLabels,:) = [];
                
                
                
                Temp = ActionTotalTensors == ActionTotalTensors(nn,3);
                ActionTensorsForgottenLabels(Temp) = repmat({NewForgottenLabels},[sum(Temp(:)),1]);
                
                TempOld = OldTrue == ActionTotalTensors(nn,3);
                if sum(TempOld(:))~=0
                    OldTrueForgottenLabels(TempOld) = repmat({NewForgottenLabels},[sum(TempOld(:)),1]);
                end
            end
            
        end
        
        if ~isempty(ActionTensorsPermute{nn,3})
            for kk = 1:numel(CTensors)
                CTensors{kk} = permute(CTensors{kk}, ActionTensorsPermute{nn,3});
            end
            CSizesUsed = CSizesUsed(:,ActionTensorsPermute{nn,3});
        end
        
        %store the output:
        
        TensorsRaw{ActionTensorsLocations(nn,3)} = CTensors;
        TensorsSizes{ActionTensorsLocations(nn,3)} = CSizesUsed;
        
    else %ActionTypeTotal(nn) == 0
        %then we are doing a kronica product
        
        
        
        %first take in what I'm going to be working with and perform the
        %moves so that the Tensor entries are in the correct order
        ATensors = TensorsRaw{ActionTensorsLocations(nn,1)};
        ASizes = TensorsSizes{ActionTensorsLocations(nn,1)};
        
        if ~isempty(ActionTensorsPermute{nn,1})
            for kk = 1:numel(ATensors)
                ATensors{kk} = permute(ATensors{kk},ActionTensorsPermute{nn,1});
            end
            ASizes = ASizes(:,ActionTensorsPermute{nn,1});
        end
        BTensors = TensorsRaw{ActionTensorsLocations(nn,2)};
        BSizes = TensorsSizes{ActionTensorsLocations(nn,2)};
        
        if ~isempty(ActionTensorsPermute{nn,2})
            for kk = 1:numel(BTensors)
                BTensors{kk} = permute(BTensors{kk},ActionTensorsPermute{nn,2});
            end
            BSizes = BSizes(:,ActionTensorsPermute{nn,2});
        end
        
        
        %delete the stuff I can forget:
        TensorsRaw(ActionTensorsLocations(nn,[DropTensor(nn,:),false])) = [];
        TensorsSizes(ActionTensorsLocations(nn,[DropTensor(nn,:),false])) = [];
        
        
        if isequal(SymmetryUsed.getSymName,'none')
            
            CSizes = [ASizes,BSizes];
            CTensors = {kron(ATensors{1},BTensors{1})};
            
            
        else
            OldForgottenLabelsA = ActionTensorsForgottenLabels{nn,1};
            TranslateA = 1:(numel(ATensors)+numel(OldForgottenLabelsA));
            Temp = TranslateA; Temp(OldForgottenLabelsA) = [];
            TranslateA(Temp) = 1:numel(ATensors);
            
            OldForgottenLabelsB = ActionTensorsForgottenLabels{nn,2};
            TranslateB = 1:(numel(BTensors)+numel(OldForgottenLabelsB));
            Temp = TranslateB; Temp(OldForgottenLabelsB) = [];
            TranslateB(Temp) = 1:numel(BTensors);
            
            List = ActionChargesList{nn};
            
            KeepA = ~any(repmat(List(:,1),[1,numel(OldForgottenLabelsA)])==repmat(OldForgottenLabelsA(:)',[size(List,1),1]),2);
            KeepB = ~any(repmat(List(:,2),[1,numel(OldForgottenLabelsB)])==repmat(OldForgottenLabelsB(:)',[size(List,1),1]),2);
            KeepC = KeepA&KeepB;
            
            CListMod = List(KeepC,:);
            [~,UseListForC,~] = unique(CListMod(:,3));
            
            CSizes = zeros([size(C_List{nn},1),size(ASizes,2)+size(BSizes,2)-2*ActionNumberContractedLegs(nn)]);
            CSizes(CListMod(UseListForC,3),:) = [ASizes(TranslateA(CListMod(UseListForC,1)),1:(end-ActionNumberContractedLegs(nn))),...
                                 BSizes(TranslateB(CListMod(UseListForC,2)),(ActionNumberContractedLegs(nn)+1):end)];
                             
            
            CTensors = cell([ActionMaxNumTensorParts(nn),1]);
            
            if ~FlagMustUsePlanar
                if ~any(List(kk,1) == OldForgottenLabelsA) && ~any(List(kk,2) == OldForgottenLabelsB)
                    for kk = 1:size(List,1)
                        CTensors{List(kk,3)} = kron(ATensors{TranslateA(List(kk,1))},BTensors{TranslateB(List(kk,2))});
                    end
                end
            else
                ListMultNumbers = ActionTensorsListMultNumbers{nn};
                if ~any(List(kk,1) == OldForgottenLabelsA) && ~any(List(kk,2) == OldForgottenLabelsB)
                    for kk = 1:size(List,1)
                        if isempty(CTensors{List(kk,3)})
                            CTensors{List(kk,3)} = kron(ATensors{TranslateA(List(kk,1))},BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                        else
                            CTensors{List(kk,3)} = CTensors{List(kk,3)} + kron(ATensors{TranslateA(List(kk,1))},BTensors{TranslateB(List(kk,2))})*ListMultNumbers(kk);
                        end
                    end
                end
            end
            
           NewForgottenLabels = zeros([1,0]);
           for ll = 1:numel(CTensors)
                if isempty(CTensors{ll})
                    NewForgottenLabels = [NewForgottenLabels,ll];
                elseif max(abs(CTensors{ll}(:)))<MinNonZero
                    NewForgottenLabels = [NewForgottenLabels,ll];
                end
            end
            
            CTensors(NewForgottenLabels) = [];
            CSizesUsed(NewForgottenLabels,:) = [];
            
            Temp = ActionTotalTensors == ActionTotalTensors(nn,3);
            ActionTensorsForgottenLabels(Temp) = repmat({NewForgottenLabels},[sum(Temp(:)),1]);
            
            TempOld = OldTrue == ActionTotalTensors(nn,3);
            if sum(TempOld(:))~=0
                OldTrueForgottenLabels(TempOld) = repmat({NewForgottenLabels},[sum(TempOld(:)),1]);
            end
        end
        
        if ~isempty(ActionTensorsPermute{nn,3})
            for kk = 1:numel(CTensors)
                CTensors{kk} = permute(CTensors{kk}, ActionTensorsPermute{nn,3});
            end
            CSizes = CSizes(:,ActionTensorsPermute{nn,3});
        end
        
        %store the output:
        
        TensorsRaw{ActionTensorsLocations(nn,3)} = CTensors;
        TensorsSizes{ActionTensorsLocations(nn,3)} = CSizes;
        
    end
    
end




































 TensorsOut = cell([numel(OldTrueLocations),1]);
for kk = 1:numel(OldTrueLocations)
    WorkingEntries = TensorsRaw{OldTrueLocations(kk)};
    
    if ~isempty(OldTruePermute{kk})
        for ll = 1:numel(WorkingEntries)
           WorkingEntries{ll} = permute(WorkingEntries{ll},OldTruePermute{kk});
        end
    end
    
    
    if ~isequal(SymmetryUsed.getSymName,'none')
        ModifiedChargeLabelsExternal = OldTrueChargesExternal{kk};
        ModifiedChargeLabelsInternal = OldTrueChargesInternal{kk};
        ModifiedMultiplicitiesInternal = OldTrueMultiplicitiesInternal{kk};
    else
        N = numel(OldTruePermute{kk});
        
        ModifiedChargeLabelsExternal = ones([1,N]);
        if ~isempty(OldTrueStructure{kk})
            ModifiedChargeLabelsInternal = ones([1,N-2]);
            ModifiedMultiplicitiesInternal = ones([1,N-2]);
        else
            ModifiedChargeLabelsInternal = ones([1,0]);
            ModifiedMultiplicitiesInternal = ones([1,0]);
        end
    end
    
    if FlagMustUsePlanar
        %now we need to play with the blocks
        List = OldTrueList{kk};
        ListMultNumbers = OldTrueListMultNumbers{kk};
        CTensors = cell([size(OldTrueChargesExternal{kk},1),1]);
        
        
        
        OldForgottenLabels = OldTrueForgottenLabels{kk};
        Translate = 1:(numel(WorkingEntries)+numel(OldForgottenLabels));
        Temp = Translate; Temp(OldForgottenLabels) = [];
        Translate(Temp) = 1:numel(WorkingEntries);
        
        for ll = 1:size(List,1)
            if ~any(List(ll,1) == OldForgottenLabels)
                if isempty(CTensors{List(ll,3)})
                    CTensors{List(ll,3)} = WorkingEntries{Translate(List(ll,1)),1}*ListMultNumbers(ll);
                else
                    CTensors{List(ll,3)} = CTensors{List(ll,3)} + WorkingEntries{Translate(List(ll,1)),1}*ListMultNumbers(ll);
                end
            end
        end
        
        NewForgottenLabels = zeros([1,0]);
        for ll = 1:numel(CTensors)
            if isempty(CTensors{ll})
                NewForgottenLabels = [NewForgottenLabels,ll];
            elseif max(abs(CTensors{ll}(:)))<MinNonZero
                NewForgottenLabels = [NewForgottenLabels,ll];
            end
        end
        
        CTensors(NewForgottenLabels) = [];
        ModifiedChargeLabelsExternal(NewForgottenLabels,:) = [];
        ModifiedChargeLabelsInternal(NewForgottenLabels,:) = [];
        ModifiedMultiplicitiesInternal(NewForgottenLabels,:) = [];
        
        WorkingEntries = CTensors;
    elseif ~isequal(SymmetryUsed.getSymName,'none')
        OldForgottenLabels = OldTrueForgottenLabels{kk};
        
        ModifiedChargeLabelsExternal(OldForgottenLabels,:) = [];
        ModifiedChargeLabelsInternal(OldForgottenLabels,:) = [];
        ModifiedMultiplicitiesInternal(OldForgottenLabels,:) = [];
    end
    
    
    TensorsOut{kk} = SymTensor.CreateRawTensors(WorkingEntries, OldTrueStructure{kk}, OldTrueChargeDirections{kk},...
                ModifiedChargeLabelsExternal, ModifiedChargeLabelsInternal, ModifiedMultiplicitiesInternal, SymmetryUsed,...
                OldTrueChargeSides{kk}(1,:),OldTrueChargeSidesInt{kk}, false,OldTrueBraiding{kk},OldTrueBraidingDirection{kk});
end

if (numel(OldTrueLocations)==1)&&FlagSingleNotCell
    TensorsOut = TensorsOut{1};
end


end

function A = ndim(B)
    A = Acore.ndim(B);
end

function bool = IsInteger(A)
    bool = abs(A-round(A))<10^-14;
end

function Bool = GTPath(In1,In2)
    NumberLength = min(size(In1,1),size(In2,1));
    
    Temp = (In1(1:NumberLength,:)'<In2(1:NumberLength,:)')-...
        (In1(1:NumberLength,:)'>In2(1:NumberLength,:)');
    
    Loc = find(Temp(:)~=0,1,'first');
    
    if isempty(Loc)
        %if they are exactly the same then this is a special case where the
        %leg from the source is being compared to the leg from one of the
        %layers it generates so we can give an exact number
        if size(In1,1)>size(In2,1)
            %1 is greater (more right) then 2
            Direction = -1;
        elseif size(In1,1)<size(In2,1)
            %1 is not greater (more right) then 2
            Direction = +1;
        else
            error('Error: The pair should not give the same path for legs')
        end
    else
        Direction = Temp(Loc);
    end
    
    if Direction == -1
        Bool = true;
    else
        Bool = false;
    end
end







function [LegOrder,FlagUserDidNotChooseLegOrder,UniqueLegsNeg, TotalTensors,TensorCopies] =...
    Broad_1_ComputingWorkingLegs(LegOrder,Legs,ChargeSideTensorDetails, ChargeDirectionsTensorDetails,StructureTensorDetails,NumberLegsTensor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STEP 1 - BROAD: Formatting LegOrder
%
%Inputs - LegOrder
%       - Legs
%       - ChargeSideTensorDetails
%       - ExternalChargeTensorDetails
%
%Outputs - LegOrder [2,4,5]
%        - FlagUserDidNotChooseLegOrder [Planar]
%        - UniqueLegsNeg [Planar]
%        - TotalTensors [2, 3, 4, 5]
%        - TensorCopies [2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlagDebug = false;



if ~iscell(Legs)
    error('Legs should be a 1xN cell of array valued entries')
end
Legs = Legs(:);

AllNumbersLegs = [];
AllNumbersSides = [];
AllNumbersDirections = [];
TensorCopies = ones([1,length(Legs)]);


for kk = 1:length(Legs)
    if ~isnumeric(Legs{kk}) || ~isreal(Legs{kk})
        error('Entry %i in Legs is not a numerical object', kk)
    end
    
    if size(Legs{kk},2)~=NumberLegsTensor(kk);
        error('Entry number %i in Legs does not line up with the corresponding tensor object', kk)
    end
    
    if any((round(Legs{kk})-Legs{kk})>10^-14 | round(Legs{kk})==0)
        error('Entry %i in Legs is not integer valued',kk)
    end
    
    Legs{kk} = round(Legs{kk});
    TensorCopies(kk) = size(Legs{kk},1);
end
TotalTensors = [0,cumsum(TensorCopies)];


    for kk = 1:length(Legs)
        ExternalCharges = zeros([1,-min(StructureTensorDetails{kk}(:))]);
        for ll = 1:-min(StructureTensorDetails{kk}(:))
            ExternalCharges(ll) = ChargeDirectionsTensorDetails{kk}(StructureTensorDetails{kk} == -ll);
        end
        
        Temp = Legs{kk}';
        TempSides = repmat(ChargeSideTensorDetails{kk}', [size(Legs{kk},1),1]);
        TempDirections = repmat(ExternalCharges, [size(Legs{kk},1),1]);
        
        AllNumbersLegs = [AllNumbersLegs, Temp(:)'];
        AllNumbersSides = [AllNumbersSides, TempSides(:)'];
        AllNumbersDirections = [AllNumbersDirections, TempDirections(:)'];
        
        ExternalChargeTensorDetails{kk} = ExternalCharges;
    end
    
    AllNumbersLegsNeg = AllNumbersLegs(AllNumbersLegs<0);
    AllNumbersSidesNeg = AllNumbersSides(AllNumbersLegs<0);
    AllNumbersDirectionsNeg = AllNumbersDirections(AllNumbersLegs<0);
    
    [AllNumbersLegsNegSorted,IndexNeg] = sort(AllNumbersLegsNeg, 'descend');
    AllNumbersSidesNegSorted = AllNumbersSidesNeg(IndexNeg);
    AllNumbersDirectionsNegSorted = AllNumbersDirectionsNeg(IndexNeg);
    
    
    if any(AllNumbersLegsNegSorted(1:(end-1)) == AllNumbersLegsNegSorted(2:end))
        error('We have double ups in the number of negative legs we are using, specifically we have extra values: %s', ...
            num2str(AllNumbersLegsNegSorted(AllNumbersLegsNegSorted(1:(end-1)) == AllNumbersLegsNegSorted(2:end))) );
    end
    
    %these will be used in Planar
    UniqueLegsNeg = AllNumbersLegsNegSorted;
    

    if ~(nargin>=4 && ~isempty(LegOrder))
        %if it is unset, input the data for contracting everything.
        LegOrder = {[],[],[]};
    end
    
%if LegOrder is numberic, then make it a cell, if it is not then make sure
%it is no larger then a 2D network, and otherwise check that it is a cell.
    if isnumeric(LegOrder)
        LegOrder = {LegOrder};
    elseif length(size(LegOrder))>2
        error('Error: LegOrder is too large');
    else
        if ~iscell(LegOrder)
            error('Error: LegOrder should be a cell or a numerical entry');
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now work out the types of entries in LegOrder (as we now know for certain
%that it is a cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TypeLegOrder = -ones(size(LegOrder)); %-1 means unchecked.
    for kk = 1:numel(LegOrder);
        if ~isempty(LegOrder{kk})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if the entry is not empty then make sure it is a 2D integer
            %array.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isnumeric(LegOrder{kk})
                error('Error: LegOrder Entries should be numerical')
            end
            if any(~IsInteger(LegOrder{kk}))
                error('Error: LegOrder Entries should be integer entries')
            end
            if length(size(LegOrder{kk}))>2||(length(size(LegOrder{kk}))>3 && size(LegOrder{kk},3) == 2)
                error('Error: LegOrder Entries should be 2D arrays (with possible detail on leg direction and charge direction on the 3rd axis)')
            end
            
            %so far this entry could be a remove, ordering or structure
            CanBeRemove = true;
            CanBeOrdering = true;
            CanBeStructure = true;
            
            %if there is a third dimension then this can't be Remove
            if size(LegOrder{kk},3)>1
                CanBeRemove = false;
            end
            
            %presence or absence of negatives means either not Remove or
            %Structure (says nothing about Ordering)
            if any(any(LegOrder{kk}<0))
                CanBeRemove = false;
            else
                CanBeStructure = false;
            end
            
            %The size in the vertical direction says something about
            %if it could be structure or ordering as Ordering must be
            %horizontal
            if size(LegOrder{kk},1) == 1
                CanBeStructure = false;
            elseif size(LegOrder{kk},1) == 2
                CanBeOrdering = false;
            else %size(LegOrder{kk},1) >2
                CanBeStructure = false;
                CanBeOrdering = false;
                if size(LegOrder{kk},2) >2
                    error('Error: LegOrder is wrong')
                end
            end
            
            
            %error catching
            if ~(CanBeOrdering || CanBeStructure||CanBeRemove)
                error('Error: Every entry in LegOrder must be either an Ordering, Structure or a Remove')
            end
            
            if CanBeStructure &&(CanBeOrdering||CanBeRemove)
                error('Affirmation Error: My Logic Failed Somehow')
            end
            TypeLegOrder(kk) = CanBeOrdering*1+CanBeRemove*2;
        else
            TypeLegOrder(kk) = 4;%set to 4 if empty
        end
        %now we have what each type of thing can be
    end
    
    FlagPossibleAmbiguity = true; %not sure which way the user input data
    FlagSuccess = false; %if we succeeded
    
    if size(LegOrder,2)>3 %if 2nd direction is too large then no ambiguity
        if size(LegOrder,1)>3 %if both too large then no solution
            error('Error: The LegOrder that was input was too large in both directions')
        end
        LegOrder = LegOrder';
        TypeLegOrder = TypeLegOrder';
        FlagPossibleAmbiguity = false;
    elseif size(LegOrder,1)>3
        %if this then we have fixed the directions
        FlagPossibleAmbiguity = false;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now do a loop to work out details
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Counter = 0;
    while Counter<=FlagPossibleAmbiguity
        
        if Counter ~= 0; %if counter isn't zero this is our second loop.
            LegOrder = LegOrder';
            TypeLegOrder = TypeLegOrder';
            FlagPossibleAmbiguity = false;
        end
        Counter = Counter+1;
        %this just loops until we have checked all combos
        
        
        %check that at most there is one of each type on each row
        
        Total0s = sum(TypeLegOrder==0,2);
        Total1s = sum(TypeLegOrder==1,2);
        Total2s = sum(TypeLegOrder==2,2);
        Total3s = sum(TypeLegOrder==3|TypeLegOrder==2|TypeLegOrder==1,2);
        Total3sOnly = sum(TypeLegOrder==3,2);
        if ~(all(Total0s<=1)&&all(Total1s<=1)&&all(Total2s<=1)&&all(Total3s<=2))
            %this orientation can't work
            if FlagDebug; disp('Error in: Totals don''t work in this orientation'); end
            continue;
        end
        
        %I know what everything must be unless there are two 3s appearing
        TempTypeLegOrder = TypeLegOrder;
        TempLegOrder = LegOrder;
        
        
        %go through each row
        for kk = 1:size(TempTypeLegOrder,1) 
            
            %if zero 3s then I know and don't have to change anything, if
            %one 3 then I want to check other details to work out which one 
            %it may be forced to be (if no ordering then it must be remove 
            %as I can't have a non-negative ordering without a remove).
            
            if Total3sOnly(kk) == 1
                if Total1s(kk) == 1
                    %then its a 2 (ordering)
                    TempTypeLegOrder(kk,TempTypeLegOrder(kk,:) == 3) = 2;
                else
                    %then it has to be a 1 (remove)
                    TempTypeLegOrder(kk,TempTypeLegOrder(kk,:) == 3) = 1;
                end
                
            %if there are two 3s then I need to work out which is which.
            elseif Total3sOnly(kk)==2
                Checking = TempLegOrder(kk,TempTypeLegOrder(kk,:) == 3);
                %see if they can be remove:
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %This section works out which of the two could be remove
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                CheckingCanBeRemove = [true,true];
                for ff = 1:2
                    Temp = Checking{ff}; 
                    %here I don't need to worry about the possibility that 
                    %it could be 3D as it could be remove (eliminating the 
                    %possibility size(Temp,3)>1
                    
                    HorizontalFine = all(Temp(1,:)<size(TotalTensors,2))&&all(Temp(1,:)>0);
                    VerticalFine = all(Temp(:,1)<size(TotalTensors,2))&&all(Temp(:,1)>0);
                   
                    HorizontalCount = true;
                    if HorizontalFine&&size(Temp,1)>1
                       HorizontalCount = all(Temp(2,:)<=(TensorCopies(Temp(1,:))));
                    end
                    VerticalCount = true;
                    if VerticalFine&&size(Temp,2)>1
                        VerticalCount = all(Temp(:,2)<=(TensorCopies(Temp(:,1))));
                    end
                    
                    if VerticalFine&&VerticalCount&&size(Temp,2)<=2
                        if size(Temp,2) == 1
                            Temp = [Temp,ones(size(Temp))];
                        end
                    elseif HorizontalFine&&HorizontalCount&&size(Temp,1)<=2
                        Temp = Temp';
                        if size(Temp,2) == 1
                            Temp = [Temp,ones(size(Temp))];
                        end
                    else
                        CheckingCanBeRemove(ff) = false;
                    end
                    
                    if ff == 1
                        Temp1 = Temp;
                    else
                        Temp2 = Temp;
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %now work out which is which
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ~(CheckingCanBeRemove(1) || CheckingCanBeRemove(2))
                    %if neither can be a remove then this orientation fails
                    if FlagDebug; disp('Error in: Neither of the options could be a remove'); end
                    continue;
                elseif (CheckingCanBeRemove(1) && CheckingCanBeRemove(2))
                    %if both could be remove
                    
                    %what are all the names for the legs (in descending order)
                    SpareLegs = AllNumbersLegsNegSorted;
                    for ll = 1:size(Temp1,1)
                        SpareLegs = [SpareLegs,Legs{Temp1(ll,1)}(Temp1(ll,2),:)];
                    end
                    SpareLegs = sort(SpareLegs,'descend');
                    %which legs are duplicates
                    SpareLegsLogic = SpareLegs(2:end) == SpareLegs(1:(end-1));
                    SpareLegsLogic = [false,SpareLegsLogic] | [SpareLegsLogic,false];
                    
                    %only keep non-duplicates (and sort in ascending order)
                    SpareLegs = sort(SpareLegs(~SpareLegsLogic),'ascend');
                    SpareLegs1 = SpareLegs;
                    
                    
                    
                    %what are all the names for the legs (in descending order)
                    SpareLegs = AllNumbersLegsNegSorted;
                    for ll = 1:size(Temp2,1)
                        SpareLegs = [SpareLegs,Legs{Temp2(ll,1)}(Temp2(ll,2),:)];
                    end
                    SpareLegs = sort(SpareLegs,'descend');
                    %which legs are duplicates
                    SpareLegsLogic = SpareLegs(2:end) == SpareLegs(1:(end-1));
                    SpareLegsLogic = [false,SpareLegsLogic] | [SpareLegsLogic,false];
                    
                    %only keep non-duplicates (and sort in ascending order)
                    SpareLegs = sort(SpareLegs(~SpareLegsLogic),'ascend');
                    SpareLegs2 = SpareLegs;
                    
                    
                    
                    
                    %check if the other entry that we are checking works
                    %with this remove (note that if they both work then we
                    %assume it was input as remove and then Ordering)
                    if isequal(sort(Checking{2}(:),'ascend'),SpareLegs1(:))
                        TempTypeLegOrder(kk,TempTypeLegOrder(kk,:) == 3) = [2,1];
                    elseif isequal(sort(Checking{1}(:),'ascend'),SpareLegs2(:))
                        TempTypeLegOrder(kk,TempTypeLegOrder(kk,:) == 3) = [1,2];
                    else
                        %neither worked so we failed
                        if FlagDebug; disp('Error in: Remove'); end
                        continue;
                    end
                    
                else %xor (CheckingCanBeRemove1 || CheckingCanBeRemove2)
                    %if only one could be remove then we can set
                    TempTypeLegOrder(kk,TempTypeLegOrder(kk,:) == 3) = 1+CheckingCanBeRemove;
                end
            end
                
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %in this case we have 0,1,2 in our listing (and 4) and so we want
        %to check that the entries actual make sense.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %base data here on TempTempWorkingOrder (as we might rotate later)
        
        FlagUserDidNotChooseLegOrder = false([size(LegOrder,1),1]);
        WorkingLegOrder = cell([size(LegOrder,1),8]);
        for kk = 1:size(LegOrder,1)
            
            PossibleChargeDirections = [];
            PossibleExternalChargeSides = [];
            
            %set the temperary working set
            WorkingSet = {zeros([0,2]),zeros([1,0]),zeros([2,0]),zeros([2,0]),zeros([1,0]),zeros([1,0]),zeros([1,0]),zeros([1,0])};
            
            %incase we jump out of the loop before fixing this up.
            WorkingLegOrder(kk,:) = WorkingSet;
            
            Loc = TempTypeLegOrder(kk,:)==2;%find the remove entry
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %check all details are correct with remove
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TempBottomNumber = -1;
            if sum(Loc)~=0
                %if there is one then
                Temp = TempLegOrder{kk,Loc};
                %if more then one direction is greater then 1
                if ((size(Temp,1)>2)&&(size(Temp,2)>2))||length(size(Temp))>2
                    %then the input was incorrect
                    if FlagDebug; disp('Error in: error in Remove - Size'); end
                    continue;
                end
                
                %check if Horizontal or vertical
                HorizontalFine = all(Temp(1,:)<size(TotalTensors,2))&&all(Temp(1,:)>0);
                VerticalFine = all(Temp(:,1)<size(TotalTensors,2))&&all(Temp(:,1)>0);
                
                HorizontalCount = true;
                if HorizontalFine&&size(Temp,1)>1
                    HorizontalCount = all(Temp(2,:)<=(TensorCopies(Temp(1,:))));
                end
                VerticalCount = true;
                if VerticalFine&&size(Temp,2)>1
                    VerticalCount = all(Temp(:,2)<=(TensorCopies(Temp(:,1))));
                end
                
                %set it up to be a row of entries and if it is only given with one row then add another
                if VerticalFine&&VerticalCount&&size(Temp,2)<=2
                    if size(Temp,2) == 1
                        Temp = [Temp,ones(size(Temp))];
                    end
                elseif HorizontalFine&&HorizontalCount&&size(Temp,1)<=2
                    Temp = Temp';
                    if size(Temp,2) == 1
                        Temp = [Temp,ones(size(Temp))];
                    end
                else
                    %then the remove is incorrect
                    if FlagDebug; disp('Error in: error in Remove - neither horizontal or vertical works'); end
                    break;
                end
                
                %store this as the first entry
                WorkingSet{1} = Temp;
            else
                Temp = WorkingSet{1};
                %because I use Temp later and WorkingSet{1} is already empty
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Here I have set WorkingSet{1}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now start working out legs (i.e. Remove)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            SpareLegs = AllNumbersLegsNegSorted;
            SpareSides = AllNumbersSidesNegSorted;
            SpareDirections = AllNumbersDirectionsNegSorted;
            for ll = 1:size(Temp,1)
                SpareLegs = [SpareLegs,Legs{Temp(ll,1)}(Temp(ll,2),:)];
                SpareSides = [SpareSides, -ChargeSideTensorDetails{WorkingSet{1}(ll,1)}]; %note that we are using the other side
                SpareDirections = [SpareDirections, -ExternalChargeTensorDetails{WorkingSet{1}(ll,1)}]; %note that we are using the other direction
            end
            
            if ~isempty(SpareLegs)
                [SpareLegs,Index] = sort(SpareLegs,'descend');
                SpareSides = SpareSides(Index);
                SpareDirections = SpareDirections(Index);
                SpareLegsLogic = SpareLegs(2:end) == SpareLegs(1:(end-1));
                SpareLegsLogic = [false,SpareLegsLogic] | [SpareLegsLogic,false];
                 
                SpareLegs = SpareLegs(~SpareLegsLogic);
                SpareSides = SpareSides(~SpareLegsLogic);
                SpareDirections = SpareDirections(~SpareLegsLogic);
                
                NumNegLegs = sum(SpareLegs<0);
                
                [SpareLegsPos, IndexPos] = sort(SpareLegs((NumNegLegs+1):end),'ascend');
                SpareLegs = [SpareLegs(1:NumNegLegs),SpareLegsPos];
                SpareSides((NumNegLegs+1):end) = SpareSides(IndexPos+NumNegLegs);
                SpareDirections((NumNegLegs+1):end) = SpareDirections(IndexPos+NumNegLegs);
            else
                %make the shape very specific
                SpareLegs = zeros([1,0]);
                SpareSides = zeros([1,0]);
                SpareDirections = zeros([1,0]);
            end
            
            %%%%%%%%%
            %check that the legs I've found are consistant with the
            %expected output
            %%%%%%%%%
            
            TempNumberBottomLegs = -1; %how many Bottom numbers
            Loc = TempTypeLegOrder(kk,:)==1;%find the (Leg) Ordering
            if sum(Loc)~=0 
                %if there is a single entry (extra entries would have been elimited already
                Temp = TempLegOrder{kk,Loc};
                
                if size(Temp,3)>1
                    PossibleExternalChargeSides = Temp(:,:,2);
                    Temp = Temp(:,:,1);
                else
                    PossibleExternalChargeSides = [];
                end
                
                if any(Temp ==0)
                    ZeroLoc = find(Temp ==0, 1,'first');
                    if numel(Temp)>ZeroLoc
                        if Temp(ZeroLoc+1)>=0
                            TempNumberBottomLegs = Temp(ZeroLoc+1);
                        end
                    end
                    Temp = Temp(1:(ZeroLoc-1));
                    if ~isempty(PossibleExternalChargeSides)
                        PossibleExternalChargeSides = PossibleExternalChargeSides(1:(ZeroLoc-1));
                    end
                end
                
                
                if ~isequal(sort(Temp(:),'ascend'),sort(SpareLegs(:), 'ascend'))
                    %then the spare legs found are inconsistant with the
                    %Leg order, ignoring actual ordering effects
                    if FlagDebug; disp('Error in: error in Remove/Ordering - they don''t agree'); end
                    break;
                end
                
                WorkingSet{2} = Temp;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Here I have WorkingSet{2}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Compute the structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Loc = TempTypeLegOrder(kk,:)==0;%find the Structure
            if sum(Loc)~=0
                Temp = TempLegOrder{kk,Loc};
                if ~isequal(sort(Temp(:),'ascend'), [(-numel(SpareLegs):-1)';(1:(numel(SpareLegs)-2))'])||numel(SpareLegs)<2||...
                    ~isequal(size(Temp),[2,numel(SpareLegs)-1])
                %if the numbers don't match structure, or if there is
                %something which is the structure and the number of legs is
                %less then 2 (and so it should be empty), or if the shape
                %is off
                    if FlagDebug; disp(['Error: Something is wrong with the structure in LegOrder at for row ',num2str(kk)]); end
                    continue;
                end
                
                if size(Temp,3)>1
                    PossibleChargeDirections = Temp(:,:,2);
                    Temp = Temp(:,:,1);
                else 
                    PossibleChargeDirections =[];
                end
                
                WorkingSet{3} = Temp;
            end
            %if there is something wrong with how the legs are connected
            %then this should be picked up when computing WorkingSet{5}
            
            if ~isempty(PossibleChargeDirections)
                PossibleExternalChargeDirections = zeros([1,size(PossibleChargeDirections,2)+1]);
                PossibleInternalChargeDirections = zeros([1,size(PossibleChargeDirections,2)-1]);
                
                for kk = 1:(size(Temp,2)+1)
                    PossibleExternalChargeDirections = PossibleChargeDirections(WorkingSet{3} == -kk);
                end
                
                for kk = 1:(size(Temp,2)-1)
                    PossibleInternalChargeDirections = PossibleChargeDirections(WorkingSet{3} == kk);
                end
            else
                PossibleExternalChargeDirections = [];
                PossibleInternalChargeDirections = [];
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Here I have WorkingSet{3}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now start working out everything else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(WorkingSet{2})
                %this might be slightly problematic when there is braiding,
                %to get around this prepare another Flag to trigger in Planar
                FlagUserDidNotChooseLegOrder(kk) = true;
                WorkingSet{2} = SpareLegs;
            else
                FlagUserDidNotChooseLegOrder(kk) = false;
            end
            
            
            
            %Now what legs are bottom legs 
            
            WorkingSet{4} = ones([2,numel(WorkingSet{2})]);
            BottomSpareSidesLegs = SpareLegs(SpareSides == -1);
            if ~isempty(BottomSpareSidesLegs)
                Temp = any(repmat(BottomSpareSidesLegs',[1,size(SpareLegs,2)])==repmat(WorkingSet{2},[size(BottomSpareSidesLegs,2),1]),1);
                WorkingSet{4}(1,Temp) = -ones([1,sum(Temp)]);
            end
            
            %now correct for legs which are flipped by the user
            if ~isempty(PossibleExternalChargeSides)
                Temp = WorkingSet{4}(1,:) ~= PossibleExternalChargeSides & (PossibleExternalChargeSides == +1 | PossibleExternalChargeSides == -1);
                WorkingSet{4}(2,:) = 1-Temp;
                WorkingSet{4}(1,:) = WorkingSet{4}(1,:).*(-1).^(Temp);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Here I have WorkingSet{4}
            %first row is the charge directions
            %second row is if they flipped
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            NegLegs = WorkingSet{4}(1,:)==-1; %bottom
            PosLegs = WorkingSet{4}(1,:)==+1; %top
            
            %this is a list of numbers 1 to max value given by
            %WorkingSet{4}
            Numbers = -(1:size(WorkingSet{4},2));
            NumNegLegs = sum(NegLegs);                        
            NumPosLegs = sum(PosLegs);
            NegLegs = Numbers(NegLegs);
            PosLegs = Numbers(PosLegs);
            
            %work out WorkingSet{3} if it wasn't already specified
            if isempty(WorkingSet{3})
                if numel(SpareLegs)<2
                    WorkingSet{3} = [];
                elseif numel(SpareLegs) ==2
                    WorkingSet{3} = [-1;-2];
                else
                    WorkingSet{3} = SymTensor.GenerateStructure([NumNegLegs,NumPosLegs]);
                end
            end
            
            
            %now that I have external sides I need to check if they
            %correspond to top and then bottom or if I needed to add additional bends to
            %make this work:
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %now work out WorkingSet{8} from WorkingSet{4} (which ones are
            %flipped)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %this element is supose to tell us how many legs are on the
            %wrong side, this combined with the sides given by
            %WorkingSet{4} tell us how many legs are on the bottom
            if TempNumberBottomLegs == -1
                WorkingSet{8} = repmat(WorkingSet{4}(1,:),[size(WorkingSet{4},2)+1,1]);
                MinimumBends = triu(ones([size(WorkingSet{4},2),size(WorkingSet{4},2)+1])')-triu(ones([size(WorkingSet{4},2),size(WorkingSet{4},2)+1]),1)';
                MinimumBends = sum(MinimumBends~=WorkingSet{8},2);
                MinimumBends = find(MinimumBends == min(MinimumBends));
                if numel(MinimumBends) > 1
                     MinimumBends = MinimumBends((MinimumBends.^2 + (numel(WorkingSet{4})+1-MinimumBends).^2)==min(MinimumBends.^2 + (numel(WorkingSet{4})+1-MinimumBends).^2));
                end
                if numel(MinimumBends) > 1
                    MinimumBends = MinimumBends(end);
                end
            else
                MinimumBends = TempNumberBottomLegs+1;
            end
            MinimumBends = [-ones([1,MinimumBends-1]),ones([1,size(WorkingSet{4},2)+1-MinimumBends])];
            WorkingSet{8} = MinimumBends~=WorkingSet{4}(1,:);
            WorkingSet{4} = WorkingSet{4}.*(-1).^repmat(WorkingSet{8},[2,1]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Here I have WorkingSet{3} and will now work out WorkingSet{5}
            %and WorkingSet{7}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %given ordering and structure we need to work out charge
            %directions and sides for internal bonds
            %note that I can assume that there are no bends now
                
            if numel(WorkingSet{2})>2
                WorkingSet{5} = ones([1,numel(WorkingSet{2})-2]);
                VertexFixed = [true, false([1,numel(WorkingSet{2})-2])];
                
                %this direction is actually the sides
                SetDirection = zeros(size(WorkingSet{3}));
                for ff = 1:-min(min(WorkingSet{3}))
                    SetDirection(WorkingSet{3} == -ff) = WorkingSet{4}(1,ff);
                end
                VertexUpdatable = all(SetDirection ~=0,1);
                
                CheckOrdering = -(1:numel(WorkingSet{2}));
                
                %now iteratively I need to update the internal bonds
                while ~all(VertexFixed)
                    
                    UpdateVertex = find(VertexUpdatable(2:end));
                    UpdateVertex = UpdateVertex(~VertexFixed(UpdateVertex+1));
                    
                    if isempty(UpdateVertex)
                        error('Affirmation Error: there should always be at least one updatable vertex');
                    end
                    
                    SignVertex = sign(sum(SetDirection(:,UpdateVertex+1),1));
                    
                    %always do the side cases last
                    if all(SignVertex == 0);
                        error('Error: Code this up')
                        %there are two of these
                        WorkingSet{5}(UpdateVertex) = -1;
                        WorkingSet{7}(UpdateVertex) = -1;
                        Numbers = 1:numel(UpdateVertex);
                        
                        for aa = Numbers(SignVertex ~=0)
                            SetDirection(UpdateVertex(aa) == WorkingSet{3}) = -1;
                            VertexFixed(UpdateVertex(aa)+1) = true;
                            
                            Loc1 = find(CheckOrdering == SetDirection(1,UpdateVertex(aa)+1));
                            Loc2 = find(CheckOrdering == SetDirection(2,UpdateVertex(aa)+1));
                            
                            if (Loc1-Loc2)==1 || (Loc1 == 1 && Loc2 == numel(CheckOrdering))|| (Loc2 == 1 && Loc1 == numel(CheckOrdering))
                                CheckOrdering(min(Loc1,Loc2)) = UpdateVertex(aa); 
                                CheckOrdering(max(Loc1,Loc2)) = [];
                            else
                                if FlagDebug; disp(['Error: the structure is not correct here']); end
                                continue;
                            end
                        end
                        
                    else
                        
                        WorkingSet{5}(UpdateVertex(SignVertex~=0)) = SignVertex(SignVertex~=0);%internal sides
                        WorkingSet{7}(UpdateVertex(SignVertex~=0)) = SignVertex(SignVertex~=0);%internal directions
                        
                        Numbers = 1:numel(SignVertex);
                        for aa = Numbers(SignVertex ~=0)
                            SetDirection(UpdateVertex(aa) == WorkingSet{3}) = SignVertex(aa);
                            VertexFixed(UpdateVertex(aa)+1) = true;
                            
                            Loc1 = find(CheckOrdering == WorkingSet{3}(1,UpdateVertex(aa)+1));
                            Loc2 = find(CheckOrdering == WorkingSet{3}(2,UpdateVertex(aa)+1));
                            
                            %check that they are neighbouring
                            if abs(Loc1-Loc2)==1
                                CheckOrdering(min(Loc1,Loc2)) = UpdateVertex(aa); 
                                CheckOrdering(max(Loc1,Loc2)) = [];
                            else
                                if FlagDebug; disp(['Error: the structure is not correct here']); end
                                continue;
                            end
                            
                        end
                        
                        VertexUpdatable = all(SetDirection ~=0,1);
                    end
                end
               
            end
            %else do nothing
            
            
            %WorkingSet{6} and WorkingSet{7}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %now work out WorkingSet{6} (like WorkingSet{4}) and fix up
            %WorkingSet{7}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            WorkingSet{6} = ones([2,numel(WorkingSet{2})]);
            BottomSpareDirectionsLegs = SpareDirections(SpareDirections == -1);
            if ~isempty(BottomSpareSidesLegs)
                Temp = any(repmat(BottomSpareSidesLegs',[1,size(SpareDirections,2)])==repmat(WorkingSet{2},[size(BottomSpareSidesLegs,2),1]),1);
                WorkingSet{6}(1,Temp) = -ones([1,sum(Temp)]);
            end
            
            if ~isempty(PossibleExternalChargeDirections)
                Temp = WorkingSet{4}(1,:) ~= PossibleExternalChargeDirections & (PossibleExternalChargeDirections == +1 | PossibleExternalChargeDirections == -1);
                WorkingSet{6}(2,:) = 1-Temp;
                WorkingSet{6}(1,:) = WorkingSet{6}(1,:).*(-1).^(Temp);
            end
            
            if ~isempty(PossibleInternalChargeDirections)
                Temp = WorkingSet{4}(1,:) ~= PossibleInternalChargeDirections & (PossibleInternalChargeDirections == +1 | PossibleInternalChargeDirections == -1);
                WorkingSet{7}(1,:) = WorkingSet{7}(1,:).*(-1).^(Temp);
            end
            
            WorkingLegOrder(kk,:) = WorkingSet;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Everything has been successfully detailed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FlagSuccess = true;
        Counter = Counter+1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %WorkingSet = {Remove, Ordering, Structure, ExternalChargeSides, InternalChargeSides, ExternalChargeDirections, InternalChargeDirections, Cup}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if FlagSuccess
        LegOrder = WorkingLegOrder;
    else
        error('Error: There is something wrong with LegOrder')
    end
end


function [UniqueLegsPos, ModifiedListIncluded,TracingTreeContractions, TracingTreeIncluded, TracingTreeTriads,...
    MaxNumberEdges, AtomicTensorTracing, AtomicTensorOutputTracing, TracedLegsContractions, NegativesTopTensors, NegativeReplacements,...
    FullContractionClosedFinalTriad, FullContractionClosedFinalContractedLegs] = ...
    Broad_2_ComputingTreeStructure(TensorCopies, ContractionOrder, TotalTensors, Legs, Crossings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STEP 2 - BROAD: Computing the Tree structure for each Contraction order
%and dropping triads and traces
%
%Inputs - TensorCopies [1]
%       - ContractionOrder [Input]
%       - TotalTensors [1]
%       - Legs [Input]
%
%Outputs - ContractionTriads [4]
%        - ContractedLegs [4]
%        - UniqueLegsPos [Planar]
%        - ModifiedListIncluded [3]
%        - ModifiedTreeTriads []
%        - ModifiedTreeTypes [3]
%        - ModifiedTreeContractions []
%        - TracingTreeIncluded [3]
%        - TracingTreeTriads [3]
%        - TracedAtomicTensors [4]
%        - MaxNumberEdges [4]
%        - AtomicTensorTracing [4]
%        - AtomicTensorOutputTracing [4]
%        - TracedLegsContractions [4]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AllNumbersLegs = [];

    for kk = 1:length(Legs)
        Temp = Legs{kk}';
        AllNumbersLegs = [AllNumbersLegs, Temp(:)'];
    end
    
    AllNumbersLegsPos = AllNumbersLegs(AllNumbersLegs>0);
    [AllNumbersLegsPosSorted,~] = sort(AllNumbersLegsPos, 'ascend');
    
    %These will be used in Planar
    UniqueLegsPos = AllNumbersLegsPosSorted(2:2:end);
    
    if nargin>=3 && ~isempty(ContractionOrder)
        if ~isnumeric(ContractionOrder) || ~isreal(ContractionOrder)
            error('The output leg order is not a numerical object')
        end
        if any((round(ContractionOrder)-ContractionOrder)>10^-14 | round(ContractionOrder)<=0)
            error('The output leg order not integer valued')
        end
        
        ContractionOrder = abs(round(ContractionOrder));
        for jj = 1:size(ContractionOrder,1)
            if ~isequal(sort(ContractionOrder(jj,:),'ascend'),UniqueLegsPos);
                error('The connection order does not agree with the number of connected legs')
            end
        end
    else
        ContractionOrder = AllNumbersLegsPosSorted(1:2:end);
    end
    
    
    for kk = 1:length(Legs)
        IndividualNumberLegsTensorPos{kk} = [0,cumsum(sum(Legs{kk}>0,2)')];
        NumberLegsTensorPos(kk) = sum(Legs{kk}(:)>0);
    end
    
    LegsAfterPos = [0,cumsum(NumberLegsTensorPos)];
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Prepare input stuff for generating contraction structure (pulling
    %apart the Leg labels and generating the DataContractions step)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LegsLabels = cell([TotalTensors(end),1]);
    for kk = 1:size(LegsLabels,1)
        LegsLabels{kk} = Legs{sum(kk>TotalTensors)}(kk-TotalTensors(sum(kk>TotalTensors)),:);
    end
    
    %there will be a number of ways to contract this network given by the
    %user here.
    
    TotalListIncluded = 0;
    
    
    IndividualTensorsNegatives = false([TotalTensors(end),1]);
    
    for kk = 1:length(Legs)
        IndividualTensorsNegatives((TotalTensors(kk)+1):TotalTensors(kk+1)) = any(Legs{kk}<0,2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for each possible contraction pathway
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    NumberInitialTraces = zeros([1,size(ContractionOrder,1)]);
    ContractionTriads = cell([1,size(ContractionOrder,1)]);
    ContractionTypes = cell([1,size(ContractionOrder,1)]);
    ContractionTriadsLegs = cell([1,size(ContractionOrder,1)]);
    ContractedLegs = cell([1,size(ContractionOrder,1)]);
    ContractionTriadsList = cell([1,size(ContractionOrder,1)]);
    FullContractionClosedFinalTriad = cell([1,size(ContractionOrder,1)]);
    FullContractionClosedFinalContractedLegs = cell([1,size(ContractionOrder,1)]);
    
    for jj = 1:size(ContractionOrder,1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %For each contraction order I need to generate a matrix giving the
        %data which is going to have no more then the number of different
        %positive legs and each entry has 7 data points:
        %    - Tensor Number 1
        %    - Which copy of that tensor
        %    - Tensor number 2
        %    - Which copy of that tensor
        %    - Unique label for input tensor 1
        %    - Unique label for input tensor 2
        %    - Are these tensors the same tensor? (true or false, i.e. 1 or 0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DataContractions = ones([7, size(UniqueLegsPos,2)]);
        
        
        %loop through each contraction
        for kk = 1:size(ContractionOrder,2)
            
            %work out the Leg name of the contraction AllNumbersLegsPos is
            %a concatenate list of legs starting from the first copy of the
            %first type of tensor, going through those copies and then onto
            %the next tensor until the last copy of the last tensor.
            kkLeg = ContractionOrder(jj,kk);
            Locations = find(AllNumbersLegsPos==kkLeg);
            if numel(Locations) ~=2
                error('problem with number of locations for value %i', kkLeg)
            end
            
            %now work out which tensor they are from LegsAfterPos indicates
            %where each different type of tensor starts.
            TensorNumber1 = sum(Locations(1)>LegsAfterPos);
            TensorNumber2 = sum(Locations(2)>LegsAfterPos);
            
            %now work out which copy they are from
            %IndividualNumberLegsTensorPos is the details of where each of
            %the copies start.
            TensorCopyNumber1 = sum((Locations(1) - LegsAfterPos(TensorNumber1))>IndividualNumberLegsTensorPos{TensorNumber1});
            TensorCopyNumber2 = sum((Locations(2) - LegsAfterPos(TensorNumber2))>IndividualNumberLegsTensorPos{TensorNumber2});
            
            
            %put in all data, the last one tells us if we are working with
            %a trace or not
            DataContractions(:,kk) = [TensorNumber1;TensorCopyNumber1;TensorNumber2; TensorCopyNumber2;...
                TensorCopyNumber1+TotalTensors(TensorNumber1);TensorCopyNumber2+TotalTensors(TensorNumber2);
                ((TensorNumber1 ~= TensorNumber2) || (TensorCopyNumber1 ~= TensorCopyNumber2))];
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now work out the tree diagram for the tensor contractions. The
        %output goal for this is:
        % - ContractionTriads
        % - ContractionType
        % - ContractionTriadsLegs
        % - ContractedLegs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %here specify top tensors that haven't been connected yet, in some
        %cases these will be kronica products rather then contractions and
        %once all contractions are taken into account then CurrentTop will 
        %list all seperate parts
        HighestTensor = [true([sum(TensorCopies),1]);false([length(UniqueLegsPos),1])];
        CurrentTop = 1:sum(TensorCopies,2);
        
        
        %prep for a loop doing all the contractions
        CurrentTensorNew = sum(TensorCopies,2)+1; %name of next output tensor
        Counter = 1;
        
        
        
        %Contraction details for later:
        % - Triad: A and B contract to give C then we get [A,B,C], trace means B = 0
        % - Type: 0 = Kronica product, 1 = trace, 2 = contraction
        % - TriadsLegs: The names of the legs of the tensor.
        % - ContractedLegs: List of Leg names in this contraction
        ContractionTriads{jj} = zeros([length(UniqueLegsPos),3]);
        ContractionTypes{jj} = zeros([length(UniqueLegsPos),1])-1;
        ContractionTriadsLegs{jj} = cell([length(UniqueLegsPos)+TotalTensors(end),1]);
        ContractionTriadsLegs{jj}(1:numel(LegsLabels)) = LegsLabels;
        ContractedLegs{jj} = cell([length(UniqueLegsPos),3]);
        
        
        
        %This is a list of the atomic tensors which are included in this
        %tensor. NumberIncluded is the number of atomic tensors
        ListIncluded{jj} = cell([sum(TensorCopies)+length(UniqueLegsPos),1]);
        for kk = 1:sum(TensorCopies);
            ListIncluded{jj}{kk} = kk;
            NumberIncluded{jj}(kk) = 1;
        end
        
        
        %first find all the traces in DataContractions
        IsATrace = DataContractions(7,:)==0;
        NumberInitialTraces(jj) = sum(IsATrace);
        
        if NumberInitialTraces(jj)>0
            TracedAtomicTensors = unique(DataContractions(5,IsATrace));
            
            for LL = TracedAtomicTensors;
                %New traces. Note: ContractionTypes can only be 1
                ContractionTypes{jj}(Counter) = 1;
                ContractionTriads{jj}(Counter, :) = [LL,0,CurrentTensorNew];
                
                
                %work out which legs are contracted (New ContractedLegs and
                %ContractionTriadsLegs entries)
                TempContractedLegs = ContractionOrder(jj,DataContractions(5,:)==LL & IsATrace);
                ContractedLegs{jj}{Counter,1} = TempContractedLegs;
                
                %work out which legs are the new contraction triad legs
                Keep1 = ~any(repmat(ContractionTriadsLegs{jj}{LL},[1,1,numel(TempContractedLegs)])==...
                                   repmat(reshape(TempContractedLegs,[1,1,numel(TempContractedLegs)]),size(ContractionTriadsLegs{jj}{LL})),3);
                ContractionTriadsLegs{jj}{CurrentTensorNew} = ContractionTriadsLegs{jj}{LL}(Keep1);
                
                HighestTensor(LL) = false;
                HighestTensor(CurrentTensorNew) = true;
                ListIncluded{jj}{CurrentTensorNew} = ListIncluded{jj}{LL};
                NumberIncluded{jj}(CurrentTensorNew) = 1;
                CurrentTop(LL) = CurrentTensorNew;
                
                %update for the next iteration:
                CurrentTensorNew = CurrentTensorNew+1;
                Counter = Counter+1;
                
            end
            
        end
        
        CompletedContractions = IsATrace;
        
        %Because I've done all the traces first I never need to worry about
        %doing a trace in an unusual direction, these are all cut out.
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now go through the contraction order and do the contractions
        %between pairs of tensors as often as possible.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kk = 1:(size(ContractionOrder,2))
            if CompletedContractions(kk)
                continue; %this contraction has already been done.
            end
            
            CurrentTensor1 = CurrentTop(DataContractions(5,kk));
            CurrentTensor2 = CurrentTop(DataContractions(6,kk));
            HighestTensor(CurrentTensor1) = false;
            HighestTensor(CurrentTensor2) = false;
            HighestTensor(CurrentTensorNew) = true;
            
            ListIncluded{jj}{CurrentTensorNew} = [ListIncluded{jj}{CurrentTensor1},ListIncluded{jj}{CurrentTensor2}];
            NumberIncluded{jj}(CurrentTensorNew) = NumberIncluded{jj}(CurrentTensor1) + NumberIncluded{jj}(CurrentTensor2);
            
            ContractionsBetweenSetsOfTensors = all(any(repmat(DataContractions(5:6,:), [1,1,numel(ListIncluded{jj}{CurrentTensorNew})])...
                       == repmat(reshape(ListIncluded{jj}{CurrentTensorNew},[1,1,numel(ListIncluded{jj}{CurrentTensorNew})]),size(DataContractions(5:6,:))),3),1);
            TempContractedLegs = ContractionOrder(jj,ContractionsBetweenSetsOfTensors & ~CompletedContractions);
            
            CompletedContractions = CompletedContractions | ContractionsBetweenSetsOfTensors;
            Keep1 = ~any(repmat(ContractionTriadsLegs{jj}{CurrentTensor1},[1,1,numel(TempContractedLegs)])==...
                               repmat(reshape(TempContractedLegs,[1,1,numel(TempContractedLegs)]),size(ContractionTriadsLegs{jj}{CurrentTensor1})),3);
            Keep2 = ~any(repmat(ContractionTriadsLegs{jj}{CurrentTensor2},[1,1,numel(TempContractedLegs)])==...
                               repmat(reshape(TempContractedLegs,[1,1,numel(TempContractedLegs)]),size(ContractionTriadsLegs{jj}{CurrentTensor2})),3);
            ContractionTriadsLegs{jj}{CurrentTensorNew} = [ContractionTriadsLegs{jj}{CurrentTensor1}(Keep1),ContractionTriadsLegs{jj}{CurrentTensor2}(Keep2)];
            
            ContractionTypes{jj}(Counter) = 2;
            ContractionTriads{jj}(Counter, :) = [CurrentTensor1,CurrentTensor2,CurrentTensorNew];
            ContractedLegs{jj}{Counter,1} = TempContractedLegs;
            
            CurrentTop(CurrentTop == CurrentTensor1) = CurrentTensorNew;
            CurrentTop(CurrentTop == CurrentTensor2) = CurrentTensorNew;
            %update for the next iteration:
            CurrentTensorNew = CurrentTensorNew+1;
            Counter = Counter+1;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now I need to worry about braiding just incase we've coupled
        %otherwise disjoint networks
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ListIncludedLegs = cell(size(ListIncluded{jj}));
        for ll = 1:numel(ListIncluded{jj})
            ListIncludedLegs{ll} = [];
             for kk = ListIncluded{jj}{ll}
                 WhichTensor = find(kk>TotalTensors,1,'last');
                 ListIncludedLegs{ll} = [ListIncludedLegs{ll}, Legs{WhichTensor}(kk-TotalTensors(WhichTensor),:)];
             end
        end
        
        HighestTensorNames = find(HighestTensor);
        SkipThisTensor = false(size(HighestTensorNames));
        AttachedTo = zeros(size(HighestTensorNames));
        
        if ~isempty(Crossings)
            for KK = 1:numel(HighestTensorNames)
                kk = HighestTensorNames(KK);
                if SkipThisTensor(KK)
                    continue;
                end
                for ll = 1:size(Crossings,1)
                    NumberInThisDisjointSet = any(ListIncludedLegs{kk} == Crossings(ll,1)) + any(ListIncludedLegs{kk} == Crossings(ll,2));
                    if NumberInThisDisjointSet == 1
                        for KK2 = (KK+1):numel(HighestTensorNames)
                            kk2 = HighestTensorNames(KK2);
                            NumberInOtherDisjointSet = any(ListIncludedLegs{kk2} == Crossings(ll,1)) + any(ListIncludedLegs{kk2} == Crossings(ll,2));
                            if NumberInOtherDisjointSet>0
                                
                                if SkipThisTensor(KK2)
                                    AttachedTo(KK) = AttachedTo(KK2);
                                    SkipThisTensor(KK) = true;
                                else
                                    AttachedTo(KK2) = kk;
                                    SkipThisTensor(KK2) = true;
                                end
                            end
                        end
                    end
                    
                    if NumberInThisDisjointSet<0 || NumberInThisDisjointSet>2
                        error('Assertion Failed: Number in this disjoint set should be betweeen 0 and 2');
                    end
                end
            end
        end
        
        %now actually contract them together:
        
        for KK = 1:numel(HighestTensorNames)
            kk = HighestTensorNames(KK);
            
            if SkipThisTensor(KK)
                continue;
            end
            
            BraidedRegions = find(AttachedTo == kk);
            
            if ~isempty(BraidedRegions)
                CurrentTensor1 = kk;
                CurrentTensor2 = HighestTensorNames(BraidedRegions(1));
                HighestTensor(CurrentTensor1) = false;
                HighestTensor(CurrentTensor2) = false;
                HighestTensor(CurrentTensorNew) = true;
                
                ContractionTypes{jj}(Counter) = 0;
                ContractionTriads{jj}(Counter, :) = [CurrentTensor1,CurrentTensor2,CurrentTensorNew];
                
                ListIncluded{jj}{CurrentTensorNew} = [ListIncluded{jj}{CurrentTensor1},ListIncluded{jj}{CurrentTensor2}];
                NumberIncluded{jj}(CurrentTensorNew) = NumberIncluded{jj}(CurrentTensor1) + NumberIncluded{jj}(CurrentTensor2);
                
                %update for the next iteration:
                CurrentTensorNew = CurrentTensorNew+1;
                Counter = Counter+1;
                for KK2 = 2:numel(BraidedRegions)
                    CurrentTensor1 = CurrentTensorNew-1;
                    CurrentTensor2 = HighestTensorNames(BraidedRegions(KK2));
                    HighestTensor(CurrentTensor1) = false;
                    HighestTensor(CurrentTensor2) = false;
                    HighestTensor(CurrentTensorNew) = true;
                    
                    ContractionTypes{jj}(Counter) = 0;
                    ContractionTriads{jj}(Counter, :) = [CurrentTensor1,CurrentTensor2,CurrentTensorNew];
                    ContractionTriadsLegs{jj}{CurrentTensorNew} = [ContractionTriadsLegs{jj}{CurrentTensor1},ContractionTriadsLegs{jj}{CurrentTensor2}];
                    
                    ListIncluded{jj}{CurrentTensorNew} = [ListIncluded{jj}{CurrentTensor1},ListIncluded{jj}{CurrentTensor2}];
                    NumberIncluded{jj}(CurrentTensorNew) = NumberIncluded{jj}(CurrentTensor1) + NumberIncluded{jj}(CurrentTensor2);
                end
            end
        end
        
        ContractionTriadsList{jj} = ListIncluded{jj}((sum(TensorCopies)+1):end);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now add the top tensors together using Kronica products.
        %Also need to construct a set of:
        %   - ModifiedTreeTriads 
        %   - ModifiedTreeTypes
        %   - ModifiedTreeContractions
        %which are the Contraction... but splitting the block into working
        %for each disjoint part of the network. For the next step we also
        %want to make sure everything is a Tree which doesn't have a fixed
        %top and want to remove the traces. These will also be modified to
        %have the same names if they turn out to be the same between
        %different contraction orders, this will speed up multiple
        %contractions by doing them both at once.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Work out the number of disjoint parts and define
        %ModifiedTreeTriads and ModifiedTreeTypes
        TopContractedTensors = find(HighestTensor);
        if jj == 1
            ModifiedTreeTriads = cell([size(ContractionOrder,1),length(TopContractedTensors)]);
            ModifiedTreeTypes = cell([size(ContractionOrder,1),length(TopContractedTensors)]);
            ModifiedTreeContractions = cell([size(ContractionOrder,1),length(TopContractedTensors)]);
            ModifiedTreeIncluded = cell([size(ContractionOrder,1),length(TopContractedTensors)]);
        end
        %Is this disjoint part a closed tensor network
        NegativesTopTensors = false([length(TopContractedTensors),1]);
        
        %Now define the Modified ContractionTriads and ContractionLegs by
        %removing all extra slots which haven't been used so far.
        ModifiedCounter = Counter-1;
        ContractionTriadsMod = ContractionTriads{jj}(1:ModifiedCounter,:);
        ContractionLegsMod = ContractedLegs{jj}(1:ModifiedCounter,:);
        
        
        %initialise the first disjoint part as the first one we computed
        OldTop = TopContractedTensors(1);
        ModifiedListIncluded = cell([1,length(TopContractedTensors)]);
        NegativeReplacements = zeros([1,length(TopContractedTensors)]);
        
        %now loop over all of the disjoint parts
        for ss = 1:length(TopContractedTensors);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Add negative numbers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %add a negative number to ModifiedListIncluded to indicate open
            %legs for this disjoint part.
            kk = (TopContractedTensors(ss));
            NegativesTopTensors(ss) = any(IndividualTensorsNegatives(ListIncluded{jj}{kk}));
            if NegativesTopTensors(ss)
                ModifiedListIncluded{ss} = [ListIncluded{jj}{kk},-ss];
            else
                ModifiedListIncluded{ss} = ListIncluded{jj}{kk};
            end
            
            %if this isn't the first tensor then we need to take the
            %Kronica product (as we need a new label for a \otimes b)
            if ss ~= 1
                %update the list from the other two 
                ListIncluded{jj}{CurrentTensorNew} = [ListIncluded{jj}{OldTop}, ListIncluded{jj}{kk}];
                NumberIncluded{jj}(CurrentTensorNew) = NumberIncluded{jj}(OldTop) + NumberIncluded{jj}(kk);
                for ll = ListIncluded{jj}{CurrentTensorNew}
                    CurrentTop(ll) = CurrentTensorNew;
                end
                
                %update OldTop so we add to the tensor product now
                OldTop = CurrentTensorNew;
                
                %update Contraction stuff (only with Kronica products so
                %only 0 for ContractionType)
                ContractionTypes{jj}(Counter) = 0;
                ContractionTriads{jj}(Counter, :) = [OldTop,kk,CurrentTensorNew];
                Counter = Counter+1;
                CurrentTensorNew = CurrentTensorNew+1;
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Seperate the Contractions into each disjoint part
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Need to construct a temporary modified one
            
            FlagContinue = true;
            Used = false([ModifiedCounter,1]);
            List = kk;
            
            %work backwards to work out what tensors are included in the
            %disjoint part
            while FlagContinue
                ListNew = [];
                for ff = List
                    
                    % if the used term is an original term then set it to true
                    Rows = ~Used&any(ContractionTriadsMod == ff,2);
                    %NOTE: because this is ContractionTraidsMod then we
                    %have already removed all the kronica products.
                    
                    Used(Rows) = true;
                    Temp = ContractionTriadsMod(Rows, ContractionTriadsMod(Rows,:)~=ff);
                    Temp = Temp(Temp~=0);%this rules out traces Counting for all cases (as zero appears multiple times)
                    ListNew = [ListNew, Temp(:)'];
                end
                
                if isempty(ListNew); FlagContinue = false; end
                List = ListNew;
                
            end
            
            
            if NegativesTopTensors(ss)
                %if there are open legs then pull apart Contractions to
                %give the following
                if numel(ModifiedListIncluded{ss}) >2
                    ModifiedTreeTriads{jj,ss} = ContractionTriadsMod(Used,:);
                    ModifiedTreeTriads{jj,ss}(ModifiedTreeTriads{jj,ss}(:,3) == kk,3) = -ss;
                    
                    NegativeReplacements(ss) = kk;
                    
                    ModifiedTreeContractions{jj,ss} = ContractionLegsMod(Used,:);
                    ModifiedTreeTypes{jj,ss} = ContractionTypes{jj}(Used);
                    ModifiedTreeIncluded{jj,ss} = ContractionTriadsList{jj}(Used);
                    
                    %update ContractionTriads as well so that we say that the
                    %final contraction is with some imaginary tensor which we
                    %always remove so that we can think of this as a closed
                    %tensor network for manipulations later
                    %ContractionTriads{jj}(ContractionTriads{jj}(:,3) == kk,3) = -ss;
                end
            else
                
                
                if numel(ModifiedListIncluded{ss}) == 1
                    %Don't actually do anything
                    
                elseif numel(ModifiedListIncluded{ss}) == 2
                    %this is an edge case where we are completely
                    %contracting two tensors together
                    
                    ModifiedTreeTriads{jj,ss} = ContractionTriadsMod(Used,:);
                    ModifiedTreeTriads{jj,ss}(ModifiedTreeTriads{jj,ss}(:,3) == kk,3) = 0;
                    ModifiedTreeContractions{jj,ss} = ContractionLegsMod(Used,:);
                    ModifiedTreeTypes{jj,ss} = ContractionTypes{jj}(Used);
                    ModifiedTreeIncluded{jj,ss} = ContractionTriadsList{jj}(Used);
                    
                    %update ContractionTriads as well so that we say that the
                    %final contraction is with some imaginary tensor which we
                    %always remove so that we can think of this as a closed
                    %tensor network for manipulations later
                    %ContractionTriads{jj}(ContractionTriads{jj}(:,3) == kk,3) = -ss;
                    
                else
                    
                    %if it was closed then we need to make the network closed
                    %for manipulation later
                    ModifiedTreeTriads{jj,ss} = ContractionTriadsMod(Used,:);
                    ModifiedTreeContractions{jj,ss} = ContractionLegsMod(Used,:);
                    
                    Temp = find(ModifiedTreeTriads{jj,ss}(:,3)==kk);
                    if numel(Temp) ~=1
                        error('Affirmation Failed: we can''t find the final tensor')
                    end
                    TagOutLegs(ss) = Temp;
                    
                    
                    MinOther = min(ModifiedTreeTriads{jj,ss}(TagOutLegs(ss),1:2));
                    MaxOther = max(ModifiedTreeTriads{jj,ss}(TagOutLegs(ss),1:2));
                    
                    ModifiedTreeTypes{jj,ss} = ContractionTypes{jj}(Used);
                    ModifiedTreeContractions{jj,ss} = ContractionLegsMod(Used,:);
                    ModifiedTreeIncluded{jj,ss} = ContractionTriadsList{jj}(Used);
                end
                
            end
        end
        
        
        %clean excess off:
        ContractionTriads{jj}((Counter):end,:) = [];
        ContractionTypes{jj}((Counter):end) = [];
        ContractedLegs{jj}((Counter):end,:) = [];
        ListIncluded{jj}((CurrentTensorNew):end) = [];
        NumberIncluded{jj}((CurrentTensorNew):end) = [];
        
        %Now we want track the total number of different tensors, this will
        %double up on the atomic tensors but will also count the different
        %levels of traces as different things
        if jj == 1
            TotalListIncluded(jj) = CurrentTensorNew-1;
        else
            TotalListIncluded(jj) = TotalListIncluded(jj-1)+CurrentTensorNew-1;
        end
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This part is to work out the contracted legs for the other two
        %contraction orders (instead of A and B goes to C we have A and C
        %or B and C).
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        for ss = 1:length(TopContractedTensors) %for every sector (every top)
            
            %this part of the code is only relevent if this sector is closed
            if ~NegativesTopTensors(ss)
                
                %we want to set the contraction that we are looking at down
                %as the second last contraction as the last contraction is
                %between an atomic tensor and the corresponding
                %environment, which does not correspond to a triad.
                
                for rr = (size(ModifiedTreeTriads{jj,ss},1)-1):-1:1
                    rrTemp = rr;
                    
                    %break if we have gotten to the traces at the begining
                    if ModifiedTreeTriads{jj,ss}(rrTemp,2) == 0
                        break;
                    end
                    
                    
                    
                    %the legs for the tensors in A, B and C (1,2 and 3)
                    LegsTemp1 = ContractionTriadsLegs{jj}{ModifiedTreeTriads{jj,ss}(rrTemp,1)};
                    LegsTemp2 = ContractionTriadsLegs{jj}{ModifiedTreeTriads{jj,ss}(rrTemp,2)};
                    LegsTemp3 = ContractionTriadsLegs{jj}{ModifiedTreeTriads{jj,ss}(rrTemp,3)};
                    
                    %Since this is a closed network then the legs which are
                    %contracted over Contractions13 Contractions23 are
                    %unions of these sets
                    
                    
                    %if any of the legs are both in the other leg we are
                    %contracting and in the path that we took
                    Contractions13 = LegsTemp1(  any(repmat(LegsTemp1,[size(LegsTemp3,2),1]) == repmat(LegsTemp3',[1,size(LegsTemp1,2)]),1)  );
                    Contractions23 = LegsTemp2(  any(repmat(LegsTemp2,[size(LegsTemp3,2),1]) == repmat(LegsTemp3',[1,size(LegsTemp2,2)]),1)  );
                    
                    %NextContractions is all the legs coming from the top
                    %tensor, this should be equal to all other legs which
                    %we sould contract at this stage.
                    if ~isequal(sort(LegsTemp3), sort([Contractions13,Contractions23]));
                        disp('  ')
                        error('Affirmation failed: I haven''t used all of NextContractions (or used too much of it)')
                    end
                    ModifiedTreeContractions{jj,ss}{rr,2} = Contractions23;%this is going right
                    ModifiedTreeContractions{jj,ss}{rr,3} = Contractions13; %this is fuseing to the bottom
                                 
                end
                
                if TagOutLegs(ss) ~= size(ModifiedTreeContractions{jj,ss},1)
                    error('Affirmation Failed: I expect the TagOutLegs parameter to agree with the calculation I have here')
                end
                
            end
            %if it is negatives then we will never need to compute these
            %contractions
        end
        
        %else there is nothing to modify as everthing is therefore a trace.
        
        
    end
    
    
    MaxNumberEdges = TotalListIncluded(end);
    
    ListIncludedAll = cell([TotalListIncluded(end),1]);
    NumberIncludedAll = zeros([TotalListIncluded(end),1]);
    ContranctionOrderTypeAll = NumberIncludedAll;
    
    
    Counter = 1;
    for jj = 1:length(ListIncluded)
        ListIncludedAll(Counter:(Counter+length(ListIncluded{jj})-1)) = ListIncluded{jj};
        NumberIncludedAll(Counter:(Counter+length(ListIncluded{jj})-1)) = NumberIncluded{jj};
        ContranctionOrderTypeAll(Counter:(Counter+length(ListIncluded{jj})-1)) = repmat(jj,[length(ListIncluded{jj}),1]);
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %We have now worked out the contraction tree, now we are going to drop 
    %triads and to do this we need to remove the traces at the beginning.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    TracingTreeTriads = cell(size(ModifiedTreeTypes));
    TracingTreeContractions = cell(size(ModifiedTreeContractions));
    TracingTreeIncluded = cell(size(ModifiedTreeIncluded));
    TracedLegsContractions = cell([1,TotalTensors(end)]);
    
    %this is just a compression trick to make it one big for loop
    Listing = [repmat(1:size(ModifiedTreeTypes,1),[1,size(ModifiedTreeTypes,2)]);...
        reshape(repmat(1:size(ModifiedTreeTypes,2),[size(ModifiedTreeTypes,1),1]),...
        [1,size(ModifiedTreeTypes,1)*size(ModifiedTreeTypes,2)])];
    
    for jjAll = Listing
        oo = jjAll(1);
        ss = jjAll(2);
        
        TracingTreeTriads{oo,ss} = ModifiedTreeTriads{oo,ss}(ModifiedTreeTypes{oo,ss}>1,:); %we can use greater then 1 because there are no kronica products here
        TracingTreeContractions{oo,ss} = ModifiedTreeContractions{oo,ss}(ModifiedTreeTypes{oo,ss}(1:(end))>1,:);
        TracingTreeIncluded{oo,ss} = ModifiedTreeIncluded{oo,ss}(ModifiedTreeTypes{oo,ss}(1:(end))>1);
        RemainingContractedLegs = ModifiedTreeContractions{oo,ss}(ModifiedTreeTypes{oo,ss}(1:(end))==1,:);
        
        
        
        %now we just relabel traces by their root name so we know what are
        %atomic triads (we can insert the corrections later at the end when we add traces).
        AtomicTensorTracing = min(ModifiedTreeTriads{oo,ss}(ModifiedTreeTypes{oo,ss}==1,[1,3]),[],2);
        AtomicTensorOutputTracing = max(ModifiedTreeTriads{oo,ss}(ModifiedTreeTypes{oo,ss}==1,[1,3]),[],2);
        
        [AtomicTensorTracing,Index] = sort(AtomicTensorTracing,'ascend');
        AtomicTensorOutputTracing = AtomicTensorOutputTracing(Index);
        RemainingContractedLegs = RemainingContractedLegs(Index);
        
        for kk = 1:length(AtomicTensorOutputTracing)
            TracingTreeTriads{oo,ss}(TracingTreeTriads{oo,ss}==AtomicTensorOutputTracing(kk)) = AtomicTensorTracing(kk);
            TracedLegsContractions{AtomicTensorTracing(kk)} = RemainingContractedLegs{kk}; %this will be the same for all contraction orderings.
        end
        
        %if it is a closed network and there are more then 2 tensors in it
        %then the top triad must be a center
        if ~NegativesTopTensors(ss) && numel(ModifiedListIncluded{ss})>2
            
            FullContractionClosedFinalTriad{ss} = TracingTreeTriads{oo,ss}(end,:);
            TracingTreeTriads{oo,ss} = TracingTreeTriads{oo,ss}(1:(end-1),:);
            FullContractionClosedFinalContractedLegs = TracingTreeContractions{oo,ss}(end,1);
            TracingTreeContractions{oo,ss} = TracingTreeContractions{oo,ss}(1:(end-1),:);
            TracingTreeIncluded{oo,ss} = ModifiedTreeIncluded{oo,ss}(ModifiedTreeTypes{oo,ss}(1:(end))>1);
        
            Loc1 = find(any(TracingTreeTriads{oo,ss} == FullContractionClosedFinalTriad{ss}(1),2));
            Loc2 = find(any(TracingTreeTriads{oo,ss} == FullContractionClosedFinalTriad{ss}(2),2));
            
            if isempty(Loc1)&&isempty(Loc2)
                error('Assertion Error: I couldn''t find either of the inputs to the top triad')
            end
            
            if ~isempty(Loc1)
                TracingTreeTraids{oo,ss} = TracingTreeTriads{oo,ss}([1:(Loc1-1),(Loc1+1):end,Loc1],:);
                TracingTreeTriads{oo,ss}(end,3) = FullContractionClosedFinalTriad{ss}(2);
                FullContractionClosedFinalTriad{ss}(1) = FullContractionClosedFinalTriad{ss}(2)+MaxNumberEdges;
                TracingTreeContractions{oo,ss} = TracingTreeContractions{oo,ss}([1:(Loc1-1),(Loc1+1):end,Loc1],:);
                TracingTreeIncluded{oo,ss} = TracingTreeIncluded{oo,ss}([1:(Loc1-1),(Loc1+1):end,Loc1],:);
            else
                TracingTreeTraids{oo,ss} = TracingTreeTriads{oo,ss}([1:(Loc2-1),(Loc2+1):end,Loc2],:);
                TracingTreeTriads{oo,ss}(end,3) = FullContractionClosedFinalTriad{ss}(1);
                FullContractionClosedFinalTriad{ss}(2) = FullContractionClosedFinalTriad{ss}(1)+MaxNumberEdges;
                TracingTreeContractions{oo,ss} = TracingTreeContractions{oo,ss}([1:(Loc2-1),(Loc2+1):end,Loc2],:);
                TracingTreeIncluded{oo,ss} = TracingTreeIncluded{oo,ss}([1:(Loc2-1),(Loc2+1):end,Loc2],:);
            end
            
        end
        
    end
end


function [NeedFullContractionClosed, EnvironmentLabelsForFullContractions, RequiredAtomicTensors, EnvironmentLabelsForAtomicTensors, NeedCenterTriadForAbove,...
    EnvironmentLabelsForCenterTriadAbove, NeedTriadsForAbove, EnvironmentLabelsForTriadsAbove, NeedTriadsForBelow, EnvironmentLabelsForTriadsBelow] = ...
    Broad_3_ComputingEnvironments(LegOrder, TotalTensors, NegativesTopTensors, TracingTreeTriads, TracingTreeIncluded, ModifiedListIncluded)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STEP 3 - BROAD: Computing the environment points which are cut in the tree
%diagrams and from this contructing the dropped tensors.
%
%Inputs - LegOrder [Input]
%       - TotalTensors [1]
%       - NegativesTopTensors [2]
%       - TracingTreeTriads [2]
%       - TracingTreeIncluded [2]
%       - ModifiedListIncluded [2]
%
%Outputs - NeedFullContractionClosed [4]
%        - EnvironmentLabelsForFullContractions [4]
%        - RequiredAtomicTensors [4]
%        - EnvironmentLabelsForAtomicTensors [4]
%        - EnvironmentLabelsForCenterTriadAbove [4]
%        - NeedTriadsForAbove [4]
%        - EnvironmentLabelsForTriadAbove [4]
%        - NeedTriadsForBelow [4]
%        - EnvironmentLabelsForTriadsBelow [4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now we have a set of LegOrder which are environments-by-2 cell, Each
    %environment can be thought of as an edge on the tree which is cut.
    %This part is to convert the environments data into a geometric
    %locations on contraction trees (note if there are disjoint regions 
    %there will have to be multiple geometric locations). 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %the first entry in the cell has to be the dropped tensors which is the
    %sets of tensors where we have to compute the removal of these tensors
    %independently 
    %
    %(Note this is NOT DroppedTriads)
    
    DroppedTensors = cell([size(LegOrder,1),1]);
    ListMinEnvironments = cell([size(LegOrder,1),size(TracingTreeTriads,2)]);
    %ListMinEnvironments has two entries, the first is the disjoint parts,
    %the second is the computations
    
    for jj = 1:size(LegOrder,1) %looping over wanted computations
        
        if ~isempty(LegOrder{jj,1})
            DroppedTensors{jj} = LegOrder{jj,1}(:,2)+TotalTensors(LegOrder{jj,1}(:,1))';
            %this sum is to work out the tensor rather then the tensor type.
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now split the environments into their seperate parts.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(LegOrder{jj,1})
            %if this is empty then all we need to remove from this disjoint
            %region is the -ss "fake" tensor
            for ss = 1:size(TracingTreeTriads,2)
                if NegativesTopTensors(ss)
                    ListMinEnvironments{jj,ss} = -ss;
                end
            end
        else    
            for ss = 1:size(TracingTreeTriads,2)
                %which dropped tensors are in this disjoint set?
                InThis = (any(repmat(DroppedTensors{jj}, [1,size(ModifiedListIncluded{ss},1)]) == repmat(ModifiedListIncluded{ss}, [size(DroppedTensors{jj},2),1]),2));
                ListMinEnvironments{jj,ss} = DroppedTensors{jj}(InThis,1);
                if NegativesTopTensors(ss);
                    %add the "fake" tensor if need be
                    ListMinEnvironments{jj,ss} = [ListMinEnvironments{jj,ss},-ss];
                end
            end
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The goal of this part of the code is to convert the inputs in the
    %TriadList and the List of Min Environments into a parameter called the
    %Node which 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    NeedFullContractionClosed = false([1,size(ListMinEnvironments,2)]); %do we need to fully contract?
    RequiredAtomicTensors = false([1,TotalTensors(end)]); %do we need to make sure to include trace of atomic (for returning)
    NeedCenterTriadForAbove = cell(size(TracingTreeTriads));
    NeedTriadsForBelow = cell(size(TracingTreeTriads));
    NeedTriadsForAbove = cell(size(TracingTreeTriads));
    
    EnvironmentLabelsForFullContractions = cell([1,size(ListMinEnvironments,2)]);
    EnvironmentLabelsForAtomicTensors = cell(size(RequiredAtomicTensors));
    EnvironmentLabelsForCenterTriadAbove = cell(size(TracingTreeTriads));
    EnvironmentLabelsForTriadsBelow = cell(size(TracingTreeTriads));
    EnvironmentLabelsForTriadsAbove = cell(size(TracingTreeTriads));
    
    for kk = 1:numel(NeedCenterTriadForAbove)
        NeedCenterTriadForAbove{kk} = false([1,3]);
        NeedTriadsForAbove{kk} = false([size(TracingTreeTriads{kk},1)-1,2]);
        NeedTriadsForBelow{kk} = false([size(TracingTreeTriads{kk},1)-1,1]);
        
        EnvironmentLabelsForCenterTriadAbove{kk} = repmat({zeros([2,0])},[1,3]);
        EnvironmentLabelsForTriadsAbove{kk} = repmat({zeros([2,0])},[size(TracingTreeTriads{kk},1)-1,2]);
        EnvironmentLabelsForTriadsBelow{kk} = repmat({zeros([2,0])},[size(TracingTreeTriads{kk},1)-1,1]);
    end
    
    for kk = 1:numel(EnvironmentLabelsForFullContractions)
        EnvironmentLabelsForFullContractions{kk} = zeros([2,0]);
    end
    
    for kk = 1:numel(EnvironmentLabelsForAtomicTensors)
        EnvironmentLabelsForAtomicTensors{kk} = zeros([2,0]);
    end
    
    
    %combine things together so there is only one for loop
    Listing = [repmat(1:size(ListMinEnvironments,1),[1,size(ListMinEnvironments,2)]);...
        reshape(repmat(1:size(ListMinEnvironments,2),[size(ListMinEnvironments,1),1]),...
        [1,size(ListMinEnvironments,1)*size(ListMinEnvironments,2)])];
    
    
    
    for jjAll = Listing
        
        
        ee = jjAll(1); %Environment
        ss = jjAll(2); %disjoint part
        Top = -1;
        FlagSucceeded = false;
        %Top:
        %    -1 - don't need anything
        %     0 - Need everything
        %     otherwise - Need bond coming out of Top tensor.
        
        for FF = ListMinEnvironments{ee,ss}
            if ~any(FF == ModifiedListIncluded{ss})
                error('Assertion Error: This element of ListMinEnvironments is not in ModifiedListIncluded');
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Check edge cases of 1 or 2 tensors in the disjoint set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if numel(ModifiedListIncluded{ss}) == 1
            %then this will be a closed network which is traced over
            if numel(ListMinEnvironments{ee,ss}) == 0
                RequiredAtomicTensors(ModifiedListIncluded{ss}) = true;
                EnvironmentLabelsForAtomicTensors{Atomic} = [EnvironmentLabelsForAtomicTensors{Atomic}, [ee;ss]];
                FlagSucceeded = true;
            end
            if numel(ListMinEnvironments{ee,ss}) > 1
                error('Assertion Error: There should have been only one ListMinEnvironments as there is only one tensor in this disjoint set')
            end
            continue;
        end
        
        if numel(ModifiedListIncluded{ss}) == 2
            if numel(ListMinEnvironments{ee,ss}) == 0
                NeedCenterTriadForAbove{1,ss}(3) = true;
                EnvironmentLabelsForCenterTriadAbove{1,ss}{3} = [EnvironmentLabelsForCenterTriadAbove{1,ss}{3},[ee;ss]];
            end
            if numel(ListMinEnvironments{ee,ss}) == 1
                Atomic = ModifiedListIncluded{ss}(ListMinEnvironments{ee,ss}~=ModifiedListIncluded{ss});
                RequiredAtomicTensors(ModifiedListIncluded{ss}) = true;
                EnvironmentLabelsForAtomicTensors{Atomic} = [EnvironmentLabelsForAtomicTensors{Atomic}, [ee;ss]];
                FlagSucceeded = true;
            end
            if numel(ListMinEnvironments{ee,ss}) > 2
                error('Assertion Error: There should have been only one ListMinEnvironments as there is only one tensor in this disjoint set')
            end
            continue;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %First check special cases, if there are zero or one entry in
        %ListMinEnvironments
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %If zero then we need to contract everything, this is only the case for closed
        %networks - CLOSED
        if isempty(ListMinEnvironments{ee,ss})
            NeedFullContractionClosed(ss) = true;
            EnvironmentLabelsForFullContractions{ss} = [EnvironmentLabelsForFullContractions{ss},[ee;ss]];
            FlagSucceeded = true;
        else
            %check that ListMinEnvironments{ee,ss} is in ModifiedListIncluded{ss}
            
            %assert
            if ~all(any(repmat(ModifiedListIncluded{ss}(:),[1,numel(ListMinEnvironments{ee,ss})]) == repmat(ListMinEnvironments{ee,ss}(:)',[numel(ModifiedListIncluded{ss}),1]),1),2)
                error('Assertion Failed: ListMinEnvironments{ee,ss} should be within ModifiedListIncluded{ss}')
            end
        end
        %NeedTopFullZeros - do we need to contract everything for this disjoint set (this is always false for open sets)
        %TopSiteEnvironmentsOrdering - For this Environment and disjoint set which contraction order will we use
        %TopSiteEnvironments - For this environment and disjoint set what is the top
        %Top - For this environment and disjoint set what tensor does the bond we are cutting come out from
        %TopSiteEnd - Do we want to keep everything below the cut (false means everything above)
        %FlagSucceeded - We have finished
        
        
        
        %if we only have one thing to remove - OPEN & CLOSED
        if numel(ListMinEnvironments{ee,ss}) == 1
            CenterTriad = size(TracingTreeTriads{1,ss},1);
            %the centertriad is the last triad
            
            if any(ListMinEnvironments{ee,ss} == TracingTreeTriads{1,ss}(CenterTriad,:))
                %if it is on the center triad then store seperately
                Loc = (ListMinEnvironments{ee,ss} == TracingTreeTriads{1,ss}(CenterTriad,:));
                NeedCenterTriadForAbove{1,ss}(Loc) = true;
                EnvironmentLabelsForCenterTriadAbove{1,ss}{Loc} = [EnvironmentLabelsForCenterTriadAbove{1,ss}{Loc},[ee;ss]];
            else
                Loc = find(any(ListMinEnvironments{ee,ss} == TracingTreeTriads{1,ss}(1:(end-1),1:2),2));
                Loc2 = ListMinEnvironments{ee,ss} == TracingTreeTriads{1,ss}(Loc,1:2);
                NeedTriadsForAbove{1,ss}(Loc) = true;
                EnvironmentLabelsForTriadsAbove{1,ss}{Loc,Loc2} = [EnvironmentLabelsForTriadsAbove{1,ss}{Loc,Loc2},[ee;ss]];
            end
            FlagSucceeded = true;
        end
        %Note that for open this will show up properly
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Next check other special cases, if there are all or all bar one 
        %entries in ListMinEnvironments
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %If we have all atomic tensors we don't actually need this disjoint 
        %region - OPEN & CLOSED
        if numel(ListMinEnvironments{ee,ss})==numel(ModifiedListIncluded{ss})
            FlagSucceeded = true;
        end
        
        
        %If we have all atomic tensors but one then we just need the atomic
        %tensor - OPEN & CLOSED
        if numel(ListMinEnvironments{ee,ss})==(numel(ModifiedListIncluded{ss})-1) && (numel(ModifiedListIncluded{ss})>2)
            %we need the second statement to make sure we don't hit any
            %edge cases of 1 or 2 tensors.
            
            %Work out what the atomic tensor is
            Atomic = ModifiedListIncluded{ss};
            for PosAtomic = Atomic(:)'
                if ~any(PosAtomic == ListMinEnvironments{ee,ss})
                    Atomic = PosAtomic;
                    break;
                end
            end
            
            RequiredAtomicTensors(Atomic) = true;
            EnvironmentLabelsForAtomicTensors{Atomic} = [EnvironmentLabelsForAtomicTensors{Atomic}, [ee;ss]];
            FlagSucceeded = true;
        end
        %Note that for open networks and we are removing a single atomic
        %tensor 
        
        
        
        %if the number of removed sites is somewhere in between, note that
        %ModifiedListIncluded has -ve legs as well
        if ~FlagSucceeded
            %check for each possible ContractOrder
            for oo = 1:size(ModifiedTreeTypes,1)
                
                %now in this tree each bond matches onto a set of atomic
                %tensors as described in ModifiedListIncluded. 
                
                
                %If we just need to tensors within this bond then we can 
                %just check if any of TracingTreeIncluded are exactly the
                %same
                if ~NegativesTopTensors(ss)
                    %for this it is worth remembering that
                    %ListMinEnvironments will be at least 2
                    for kk = 1:(numel(TracingTreeIncluded{oo,ss}))
                        SetOutputs = sort(TracingTreeIncluded{oo,ss}{kk}(:));
                        if isequal(SetOutputs, sort(ListMinEnvironments{ee,ss}(:)))
                            
                            %now need to find the other tensor source this
                            %is exactly the same as for an atomic tensor
                            
                            ImportantPAtomicTensor = TracingTreeTriads{oo,ss}(kk,3);
                            
                            CenterTriad = size(TracingTreeTriads{oo,ss},1);
                            %the centertriad is the last triad
                            
                            if any(ImportantPAtomicTensor == TracingTreeTriads{oo,ss}(CenterTriad,:))
                                %if it is on the center triad then store seperately
                                Loc = (ImportantPAtomicTensor == TracingTreeTriads{oo,ss}(CenterTriad,:));
                                NeedCenterTriadForAbove{oo,ss}(Loc) = true;
                                EnvironmentLabelsForCenterTriadAbove{oo,ss}{Loc} = [EnvironmentLabelsForCenterTriadAbove{oo,ss}{Loc},[ee;ss]];
                            else
                                Loc = find(any(ImportantPAtomicTensor == TracingTreeTriads{1,ss}(1:(end-1),1:2),2));
                                Loc2 = ImportantPAtomicTensor == TracingTreeTriads{oo,ss}(Loc,1:2);
                                NeedTriadsForAbove{oo,ss}(Loc) = true;
                                EnvironmentLabelsForTriadsAbove{oo,ss}{Loc,Loc2} = [EnvironmentLabelsForTriadsAbove{oo,ss}{Loc,Loc2},[ee;ss]];
                            end
                            
                            FlagSucceeded = true;
                        end
                    end
                end
                
                
                %If we just need to tensors within this bond then we can 
                %just check if any of TracingTreeIncluded are exactly the
                %opposite
                if ~FlagSucceeded
                    for kk = 1:(numel(TracingTreeIncluded{oo,ss})-1)
                        
                        %Set everything excluding TracingTreeIncluded
                        SetOutputs = sort(ModifiedListIncluded{ss}(:));
                        for nn = numel(SetOutputs):-1:1
                            if any(TracingTreeIncluded{oo,ss}{kk} == SetOutputs(nn))
                                SetOutputs(nn) = [];
                            end
                        end
                        
                        %now check
                        if isequal(SetOutputs, sort(ListMinEnvironments{ee,ss}(:)))
                            NeedTriadsForBelow{oo,ss}(kk) = true;
                            EnvironmentLabelsForTriadsBelow{oo,ss}{kk} = [EnvironmentLabelsForTriadsBelow{oo,ss},[ee;ss]];
                            FlagSucceeded = true;
                        end
                        
                    end
                end
                
                %Now if I have succeeded at finding a solution for this
                %environment within this disjoint set then jump out of this
                %loop
                if FlagSucceeded
                    break;
                end
                
            end
        end
        
        %if I never found a successful case then throw an error
        if ~FlagSucceeded
            error('We didn''t find the correct contraction for section %i in environment %i',ss,ee)
        end
        
    end
    
end


function [OldTrue, OldTrueLocations, ActionTotalTensors, ActionTensorsLocations, ActualLegsTotal, DropTensor, ActionTypeTotal, LegsContractedTotal, EnvironmentLegs] = ...
    Broad_4_ComputingContractionOrder(NeedFullContractionClosed, EnvironmentLabelsForFullContractions, RequiredAtomicTensors, EnvironmentLabelsForAtomicTensors, ...
    NeedCenterTriadForAbove, EnvironmentLabelsForCenterTriadAbove, NeedTriadsForAbove, EnvironmentLabelsForTriadsAbove, NeedTriadsForBelow, EnvironmentLabelsForTriadsBelow, ...
    LegOrder, TotalTensors, TracingTreeTriads, TracingTreeContractions, MaxNumberEdges, ModifiedListIncluded, AtomicTensorTracing, AtomicTensorOutputTracing,...
    TracedLegsContractions, Legs, NegativeReplacements,FullContractionClosedFinalTriad, FullContractionClosedFinalContractedLegs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STEP 4 - BROAD: Compute the order of the contractions when optimising for
%the fewest number of operations when we need to compute particular
%environments
%
%Inputs - NeedFullContractionClosed [3]
%       - EnvironmentLabelsForFullContractions [3]
%       - RequiredAtomicTensor [3]
%       - EnvironmentLabelsForAtomicTensors [3]
%       - NeedCenterTriadForAbove [3]
%       - EnvironmentLabelsForCenterTriadAbove [3]
%       - NeedTriadsForAbove [3]
%       - EnvironmentLabelsForTriadsAbove [3]
%       - NeedTriadsForBelow [3]
%       - EnvironmentLabelsForTriadsBelow [3]
%       - LegOrder [Input]
%       - TotalTensors [1]
%       - TracingTreeTriads [2]
%       - TracingTreeContractions [2]
%       - MaxNumberEdges [2]
%       - ModifiedListIncluded [2]
%
%Outputs -
%        - OldTrue
%        - OldTrueLocations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This first part works out what the contraction order has to be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    ActionTypeTemp = cell(size(TracingTreeTriads));
    ActionTensorsTemp = cell(size(TracingTreeTriads));
    EnvironmentLabelsOutputTemp = cell(size(TracingTreeTriads));
    LegsContractedTemp = cell(size(TracingTreeTriads));
    
    Listing = [repmat(1:size(TracingTreeTriads,1),[1,size(TracingTreeTriads,2)]);...
        reshape(repmat(1:size(TracingTreeTriads,2),[size(TracingTreeTriads,1),1]),...
        [1,size(TracingTreeTriads,1)*size(TracingTreeTriads,2)])];
    
    
    for jjAll = Listing
        
        oo = Listing(1); %Contraction ordering
        ss = Listing(2); %Disjoint set number
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %The key idea here is to drop triads from the tree and from this we
        %can redefine contracted sets of atomic tensors as psuedo atomic,
        %we drop triads until the pseudo atomic tensors are barred (have an
        %environment corresponding to a cut just above it).
        %
        %By doing this we work out what contractions can always be combined
        %together.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        PAtomicTensors = ModifiedListIncluded{ss}; %A list of what atomic tensors there are (we will make psudeo atomic tensors out of this by dropping triads).
        PAtomicTensors = PAtomicTensors(PAtomicTensors>0); %Ignore the "fake" tensor if the disjoint set is open
        PAtomicTensorsUsed = false(size(PAtomicTensors)); %A list of atomic tensors edges.
        
        
        TriadList = TracingTreeTriads{oo,ss}(1:(end-1),:); %the list of triads
        CenterTriadDetails = TracingTreeTriads{oo,ss}(end,:);
        ContractionList = TracingTreeContractions{oo,ss}; %the list of contracted legs associated to each triad
        IsDropped = false([size(TriadList,1),1]); %boolean: Have I dropped this triad so far?
        
        
        %ASSERTION
        if sum(unique(TriadList(:))<0)>1
            error('Assertion Failed: There is the wrong number of negative numbers');
        end
        
        %Now start the loop
        FlagContinue = true;
        
        TempTriadBlocks = NeedTriadsForAbove{oo,ss};
        ListExtraBlocks = [];
        
        while FlagContinue
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This loop is going through all psuedo-atomic tensors that we
            %have defined so far.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ll = 0;
            while ll < length(PAtomicTensors)
                
                %find an unused psuedo atomic tensor
                ll = ll+1;
                kk = PAtomicTensors(ll);
                
                if PAtomicTensorsUsed(ll)
                    continue;
                end
                
                if any(CenterTriadDetails == kk)
                    PAtomicTensorsUsed(ll) = true;
                    continue;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Now that we have chosen one find which triad this is
                %located in which isn't already dropped
                %
                % - Important to note here that we started with the
                % ModifiedListIncluded for this disjoint set which is the
                % boundary of this tree, if we drop a triad we need to add
                % a new tensor to this set of psuedo-Atomic tensors as the
                % boundary moves inward
                % - The boundary stops moving in if we hit something which
                % is barred where we need to keep what is above. We don't
                % need to worry about what is below because the fact we are
                % dropping the triad fixes a direction which will include
                % this cut as something that will be computed.
                % - There is an edge case where which is what happens at
                % the last point on the final triad for a closed network
                % (this doesn't effect an open network)
                % - There is an edge case for what to do if we have only
                % one bar
                % - The final triad should be treated as impassible by the
                % boundary, in the case of the open network we don't need
                % to contract, in the case of the closed network this will
                % retract to a spine if needed.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                TriadNumber = find(any(TriadList(:,1:2) == kk,2)&~IsDropped);
                %ASSERTED
                if length(TriadNumber)>1
                    error('TriadNumber is too long')
                end
                
                Bottom1 = kk;
                
                %check if the other Two in Triad are Bottom
                
                OtherTriad = TriadList(TriadNumber, :);
                
                BottomLoc = find(TriadList(TriadNumber, 1:2)==kk);
                OtherTriadLabel1 = OtherTriad(1,3-BottomLoc);
                
                
                if ~TempTriadBlocks(TriadNumber) && any(OtherTriadLabel1==PAtomicTensors)
                    %if we can drop this traid
                    PAtomicTensors = [PAtomicTensors,OtherTriad(1,3)];
                    PAtomicTensorsUsed = [PAtomicTensorsUsed, false];
                    PAtomicTensorsUsed((PAtomicTensors == kk) | (PAtomicTensors == OtherTriadLabel1)) = true;
                    IsDropped(TriadNumber) = true;
                    break;
                elseif TempTriadBlocks(TriadNumber)
                    %if we can't drop this triad
                    PAtomicTensorsUsed((PAtomicTensors == kk) | (PAtomicTensors == OtherTriadLabel1)) = true;
                    ListExtraBlocks = [ListExtraBlocks, OtherTriad(1,3)];
                elseif any(OtherTriadLabel1 == ListExtraBlocks)
                    %then we also can't drop this triad
                    PAtomicTensorsUsed(PAtomicTensors == kk) = true;
                    ListExtraBlocks = [ListExtraBlocks, OtherTriad(1,3)];
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if there are no more possible update then exit this loop, this
            %is guarenteed to converge as each loop will remove one
            %possible element from this list. [check this]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if all(PAtomicTensorsUsed)
                FlagContinue = false;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Here we are listing which environments are associated to which
        %triad outputs. It is important to learn what environments are
        %now associated to the psuedo atomic tensors or are hidden in the
        %dropped triads. 
        %
        %The left over triads form what is called spines. It is also
        %important to work out what environments are still part of the
        %spines and where they are now.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        DroppedTriads = TriadList(IsDropped,:);
        SpineTriads = TriadList(~IsDropped,:);
        
        ForAction_NeedTraidsForAbove = NeedTriadsForAbove{oo,ss}(~IsDropped,:);
        ForAction_NeedTraidsForBelow = NeedTriadsForBelow{oo,ss}(~IsDropped,:);
        ForAction_EnvironmentLabelsForTraidsForAbove = EnvironmentLabelsForTriadsAbove{oo,ss}(~IsDropped,:);
        ForAction_EnvironmentLabelsForTraidsForBelow = EnvironmentLabelsForTriadsBelow{oo,ss}(~IsDropped,:);
        ForAction_LegContractions = TracingTreeContractions{oo,ss}([~IsDropped;false],:);
        
        ForAction_NeedCenterTriadForAbove = NeedCenterTriadForAbove{oo,ss};
        ForAction_EnvironmentLabelsForCenterTriadAbove = EnvironmentLabelsForCenterTriadAbove{oo,ss};
        ForAction_CenterLegContractions = TracingTreeContractions{oo,ss}(end,:);
        
        ForAction_DroppedTriadNeededForBelow = NeedTriadsForBelow{oo,ss}(IsDropped,:);
        ForAction_DroppedTriadLegContractions = TracingTreeContractions{oo,ss}(IsDropped);
        ForAction_EnvironmentLabelsForDroppedTriads = EnvironmentLabelsForTriadsBelow{oo,ss}([IsDropped;false],:);
        
        if oo == 1
            ForAction_NeedFullContractionClosed = NeedFullContractionClosed(ss);
            ForAction_EnvironmentLabelsForFullContractionClosed = EnvironmentLabelsForFullContractions{ss};
        else
            ForAction_NeedFullContractionClosed = false;
            ForAction_EnvironmentLabelsForFullContractionClosed = [];
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now that we have removed all the triads we needed to drop we can
        %reformulate the triads in terms of spines which are columns of
        %triads where the right input triad is a psuedo-atomic tensor and
        %the left input is the previous triad's output. In the case of the
        %bottom triad of the column the left input is just another
        %pseudo-atomic tensor.
        %
        %The primary spine has the top triad being associated to the centre
        %triad as defined in the dropped triads section. Spines may come
        %off the primary spine as the psuedo-atomic tensors of the right
        %input triads (though not the left input of the bottom triad as
        %that would just be part of the original spine). Any spines which
        %are secondary come off other spines (both primary and secondary)
        %and have an open triad output at the top.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %This approach is a simplistic approach designed to make sure that
        %it can easily be decoded when read, here we don't have to worry
        %about the time or memory cost, the rest of the SymCon function is
        %likely to make these costs negligable.
        
        %make a list of what the sources are
        CentersValid = false([1,3]);
        PAtomicTensors = [PAtomicTensors, ModifiedListIncluded{ss}(ModifiedListIncluded{ss}<0)];%add all negatives
        
        
        CentersValid(1) = ~any(CenterTriadDetails(1) == PAtomicTensors);
        CentersValid(2) = ~any(CenterTriadDetails(2) == PAtomicTensors);
        CentersValid(3) = ~any(CenterTriadDetails(3) == PAtomicTensors);
        
        Sources = CenterTriadDetails(CentersValid);
        SplitsFrom = zeros([3,numel(Sources)]);
        
        if numel(Sources)~=0
            SplitsFrom(3,:) = find(CentersValid);
        end
        
        SpineTriadLength = zeros([1,0]);
        SpineTriadDetails = cell([1,0]);
        SpineSwapDetails = cell([1,0]);
        
        
        Counter = 0;
        while Counter<numel(Sources)
            
            Counter = Counter+1;
            CurrentSourceOrig = Sources(Counter);
            
            %now starting from the source work out what triads are in the
            %list
            NewAdditionalSources = [];
            CurrentSource = Sources(Counter);
            OffShootValid = false([1,2]);
            CurrentTriadList = [];
            CurrentSwapTriad = [];
            
            
            
            Counter2 = 0;
            Branches = [];
            FlagContinue = true;
            while FlagContinue
                Counter2 = Counter2+1;
                %Find the triad of the current source
                CurrentTriadNumber = find(SpineTriads(:,3) == CurrentSource);
                CurrentOffShoots = SpineTriads(CurrentTriadNumber,1:2);
                CurrentTriadList = [CurrentTriadList, CurrentTriadNumber];
                CurrentSwapTriad = [CurrentSwapTriad, false];
                
                %Work out if the triad terminates here, otherwise make one
                %of the non-psuedo atomic tensors the next current source
                %and update
                OffShootValid(1) = ~any(CurrentOffShoots(1) == PAtomicTensors);
                OffShootValid(2) = ~any(CurrentOffShoots(2) == PAtomicTensors);
                
                if all(OffShootValid)
                    NewAdditionalSources = [NewAdditionalSources, CurrentOffShoots(2)];
                    CurrentSource = CurrentOffShoots(1);
                    
                    Branches = [Branches,Counter2];
                    
                elseif any(OffShootValid)
                    
                    if OffShootValid(1)
                        CurrentSource = CurrentOffShoots(1);
                    else
                        CurrentSource = CurrentOffShoots(2);
                        CurrentSwapTriad(end) = true;
                    end
                else
                    FlagContinue = false;
                end
            end
            
            %what are the 
           SpineTriadLength = [SpineTriadLength, numel(CurrentTriadList)];
           SpineTriadDetails{Counter} = CurrentTriadList(1,end:-1:1);
           SpineSwapDetails{Counter} = CurrentSwapTriad(1,end:-1:1);
           Sources = [Sources, NewAdditionalSources];
           
           if numel(NewAdditionalSources) > 0
               SplitsFrom = [SplitsFrom,  [repmat([CurrentSourceOrig; Counter], [1,numel(NewAdditionalSources)]); (SpineTriadLength(end)+1-Branches)] ];
           end
           
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Now to characterise the spines we have generated:
        % - SplitsFrom
        % - SpineTriadDetails
        % - SpineSwapDetails
        % - Sources
        %
        %
        %Characterised by:
        % - SplitsFrom: a 3-by-N array of entries where the first row
        % indicates which triad this spine splits from, the second 
        % indicates which spine number this spine is splitting from 
        % (0 indicates the center), finally the third row indicates 
        % which triad this was from start with 1 at the top triad. If the 
        % spine splits from the center then we use 0 to indicate the spine 
        % and the element number is which entry off the center triad in 
        % entries 1-3.
        %
        % - SpineTriadDetails: A cell where the numbers indicating which 
        % rows of TriadList will makes up this spine is stored for each 
        % spine, note that it will always be from bottom to top of the 
        % spine as that is the order that the generated contraction order 
        % will HAVE to be in.
        %
        % - SpineTriadLength: The length of the entries in
        % SpineTriadDetails (As an array)
        %
        % - SpineSwapDetails: A cell where each element in the list
        % indicates if the triad was swapped compared to the list in
        % SpineTriadDetails.
        %
        % - Sources: The names of the tensors that are joining the spine to
        % its source spine/center triad.
        %
        % for later logic it is important to note that
        % SpineTriadDetails(2,:) is highly useful because if the entry kk
        % is mm then the spine came off the NN-mm+1 triad
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now need to convert to ActionTensors, Action and Keep.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %first take into account all the Spines that need to be 
        
        ContractUpTo = zeros(size(SpineTriadLength));
        ContractDownTo = SpineTriadLength+1;
        
        ModifiedTriadsForAbove = cell(size(SpineTriadDetails));
        ModifiedEnvironmentLabelsForAbove = cell(size(SpineTriadDetails));
        ModifiedEnvironmentLabelsForBelow = cell(size(SpineTriadDetails));
        ModifiedLegContractions = cell(size(SpineTriadDetails));
        ModifiedTriadList = cell(size(SpineTriadDetails));
        ContractOff = cell(size(SpineTriadDetails));
        
        for kk = 1:numel(SpineTriadDetails);
            ModifiedTriadList{kk} = TriadList(SpineTriadDetails{kk},:);
            
            ModifiedLegContractions{kk} = ForAction_LegContractions(SpineTriadDetails{kk},:);
            ModifiedTriadsForAbove{kk} = ForAction_NeedTraidsForAbove(SpineTriadDetails{kk},:);
            ModifiedEnvironmentLabelsForAbove{kk} = ForAction_EnvironmentLabelsForTraidsForAbove(SpineTriadDetails{kk},:);
            ModifiedEnvironmentLabelsForBelow{kk} = ForAction_EnvironmentLabelsForTraidsForBelow(SpineTriadDetails{kk},:);
            
            ContractOff{kk} = false([1,SpineTriadLength(kk)]);
            
            if sum(SpineSwapDetails{kk})>0
                ModifiedTriadList{kk}(SpineSwapDetails{kk},1:2) = ModifiedTriadList{kk}(SpineSwapDetails{kk},[2,1]);
                ModifiedTriadsForAbove{kk}(SpineSwapDetails{kk}, 1:2) = ModifiedTriadsForAbove{kk}(SpineSwapDetails{kk}, [2,1]);
                ModifiedEnvironmentLabelsForAbove{kk}(SpineSwapDetails{kk},1:2) = ModifiedEnvironmentLabelsForAbove{kk}(SpineSwapDetails{kk},[2,1]);
                ModifiedLegContractions{kk}(SpineSwapDetails{kk},2:3) = ModifiedLegContractions{kk}(SpineSwapDetails{kk},[3,2]);
            end
            
            if sum(ModifiedTriadsForAbove{kk}(:,1))>0
                ContractDownTo(kk) = find(ModifiedTriadsForAbove{kk}(:,1),1,'first');
            end
            
            if sum(ForAction_NeedTraidsForBelow(SpineTriadDetails{kk}))>0
                ContractUpTo(kk) = find(ForAction_NeedTraidsForBelow(SpineTriadDetails{kk}),1,'last');
            end
            
            if sum(ModifiedTriadsForAbove{kk}(:,2))>0
                ContractOff{kk} = ModifiedTriadsForAbove{kk}(:,2);
                ContractDownTo(kk) = min(ContractDownTo(kk),find(ModifiedTriadsForAbove{kk}(:,2),1,'first')+1);
                ContractUpTo(kk) = max(ContractUpTo(kk),find(ModifiedTriadsForAbove{kk}(:,2),1,'last')-1);
            end
            
        end
        
        %first work out the ordering of the spines, SplitsFrom tells us
        %this
        
        SpineOrder = [];
        GetSpinesOut = 0;
        
        while ~isempty(GetSpinesOut)
            NextSpines = find(SplitsFrom(2,:) == GetSpinesOut(1));
            SpineOrder = [SpineOrder, NextSpines];
            GetSpinesOut = [NextSpines, GetSpinesOut(2:end)];
        end
        %SpineOrder = SpineOrder(end:-1:1);
        
        %work out the contract down going from bottom to top
        for kk = SpineOrder(end:-1:1);
            if ContractDownTo(kk) ~= SpineTriadLength(kk)+1
                %then the one above it must be contracted down to
                if SplitsFrom(2,kk) == 0
                    %then it is a center case
                    ForAction_NeedCenterTriadForAbove(SplitsFrom(3,kk)) = true;
                else
                    AboveSpine = SplitsFrom(2,kk);
                    ContractOff{AboveSpine}(SplitsFrom(3,kk)) = true;
                    
                    ContractDownTo(AboveSpine) = min(ContractDownTo(AboveSpine),SplitsFrom(3,kk)+1);
                    ContractUpTo(AboveSpine) = max(ContractUpTo(AboveSpine),SplitsFrom(3,kk)-1);
                end
                
            end
        end
        
        %add ContractUp from the center
        
        NeededFromCenter = false([1,3]);
        if ForAction_NeedCenterTriadForAbove(1)
            NeededFromCenter(2) = true;
            NeededFromCenter(3) = true;
        end
        if ForAction_NeedCenterTriadForAbove(2)
            NeededFromCenter(1) = true;
            NeededFromCenter(3) = true;
        end
        
        if ForAction_NeedFullContractionClosed
            ForAction_NeedCenterTriadForAbove(3) = true;
            NeededFromCenter(3) = true;
        end
        
        if ForAction_NeedCenterTriadForAbove(3)
            NeededFromCenter(1) = true;
            NeededFromCenter(2) = true;
        end
        
        
        if NeededFromCenter(1)
            Loc = find(SplitsFrom(2,:) == 0 & SplitsFrom(3,:) == 1);
            if ~isempty(Loc)
                ContractUpTo(Loc) = SpineTriadLength(Loc);
            end
        end
        if NeededFromCenter(2)
            Loc = find(SplitsFrom(2,:) == 0 & SplitsFrom(3,:) == 2);
            if ~isempty(Loc)
                ContractUpTo(Loc) = SpineTriadLength(Loc);
            end
        end
        if NeededFromCenter(3)
            Loc = find(SplitsFrom(2,:) == 0 & SplitsFrom(3,:) == 3);
            if ~isempty(Loc)
                ContractUpTo(Loc) = SpineTriadLength(Loc);
            end
        end
        
        
        %now propegate from top to bottom:
        for kk = SpineOrder
            NeededSpines = find(SplitsFrom(2,:) == kk &( SplitsFrom(3,:)>=ContractDownTo(kk) | SplitsFrom(3,:)<=ContractUpTo(kk) ));
            if ~isempty(NeededSpines)
                for kk2 = NeededSpines
                    ContractUpTo(kk2) = SpineTriadLength(kk2);
                end
            end
        end
        
        
        
        
        
        
        
        %Now I know exactly what things need to be contracted and I can
        %just order it fully
        
        ActionTensorsTemp{oo,ss} = zeros([0,3]);
        EnvironmentLabelsOutputTemp{oo,ss} = cell([0,1]);
        LegsContractedTemp{oo,ss} = cell([0,1]);
        
        %ActionTensorsTempOrdering{oo,ss} = zeros([0,2]);
        
        
        %first Contract all going up starting from the bottom spines
        
        for kk = SpineOrder(end:-1:1)
            ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; ModifiedTriadList{kk}(1:ContractUpTo(kk),:)];
            EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}; ModifiedEnvironmentLabelsForBelow{kk}(1:ContractUpTo(kk),:)];
            LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; ModifiedLegContractions{kk}(1:ContractUpTo(kk),1)];
        end
        
        
        
        
        
        %add the contractions with the center
        
        if ForAction_NeedCenterTriadForAbove(1)
            ActionTensorsTempTemp = CenterTriadDetails(1,[2,3,1]);
            if any(ActionTensorsTempTemp < 0)
                ActionTensorsTempTemp(ActionTensorsTempTemp<0) = NegativeReplacements(ss);
            end
            ActionTensorsTempTemp = ActionTensorsTempTemp+[0,0,1]*MaxNumberEdges;
            
            ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; ActionTensorsTempTemp];
            EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}; ForAction_EnvironmentLabelsForCenterTriadAbove(1)];
            LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; ForAction_CenterLegContractions(1,2)];
        end
        
        if ForAction_NeedCenterTriadForAbove(2)
            ActionTensorsTempTemp = CenterTriadDetails(1,[1,3,2]);
            if any(ActionTensorsTempTemp < 0)
                ActionTensorsTempTemp(ActionTensorsTempTemp<0) = NegativeReplacements(ss);
            end
            ActionTensorsTempTemp = ActionTensorsTempTemp+[0,0,1]*MaxNumberEdges;
            
            ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; ActionTensorsTempTemp];
            EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}; ForAction_EnvironmentLabelsForCenterTriadAbove(2)];
            LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; ForAction_CenterLegContractions(1,3)];
        end
        
        if ForAction_NeedCenterTriadForAbove(3)
            ActionTensorsTempTemp = CenterTriadDetails(1,[1,2,3]);
            if any(ActionTensorsTempTemp < 0)
                ActionTensorsTempTemp(ActionTensorsTempTemp<0) = NegativeReplacements(ss);
            end
            ActionTensorsTempTemp = ActionTensorsTempTemp+[0,0,1]*MaxNumberEdges;
            
            ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; ActionTensorsTempTemp];
            EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}; ForAction_EnvironmentLabelsForCenterTriadAbove(3)];
            LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; ForAction_CenterLegContractions(1,1)];
        end
        
        if ForAction_NeedFullContractionClosed
            ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; FullContractionClosedFinalTriad{ss}];
            EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}; ForAction_EnvironmentLabelsForFullContractionClosed];
            LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; FullContractionClosedFinalContractedLegs{ss}];
        end
        
        
        
        
        
        
        %now contract going down from the top
        
        for kk = SpineOrder
            ActionTensorsTempTemp = ModifiedTriadList{kk}(ContractOff{kk},[1,3,2]);
            ActionTensorsTempTemp = [ModifiedTriadList{kk}(end:-1:ContractDownTo(kk),[2,3,1]); ActionTensorsTempTemp(end:-1:1,:)];
            
            EnvironmentLabelsOutputTempTemp = ModifiedEnvironmentLabelsForAbove{kk}(ContractOff{kk},2);
            EnvironmentLabelsOutputTempTemp = [ModifiedEnvironmentLabelsForAbove{kk}(end:-1:ContractDownTo(kk),1);EnvironmentLabelsOutputTempTemp(end:-1:1,1)];
            
            LegsContractedTempTemp = ModifiedLegContractions{kk}(ContractOff{kk},3);
            LegsContractedTempTemp = [ModifiedLegContractions{kk}(end:-1:ContractDownTo(kk),2);LegsContractedTempTemp(end:-1:1,1)];
            
            %adjust for taking the edge in the opposite direction
            ActionTensorsTempTemp = ActionTensorsTempTemp + repmat(MaxNumberEdges*[0,1,1], [size(ActionTensorsTempTemp,1),1]);
            
            InsertedMidContract = SpineTriadLength(kk)-find(ContractOff{kk})+1;
            CorrectOrdering = 1:(SpineTriadLength(kk)-ContractDownTo(kk)+1);
            
            if ~isempty(InsertedMidContract)
                Counter = numel(CorrectOrdering)+1;
                for kk2 = InsertedMidContract
                    CorrectOrdering = [CorrectOrdering(1:(kk2-1)),Counter, CorrectOrdering(kk2:end)];
                    Counter = Counter+1;
                end
            end
            
            ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; ActionTensorsTempTemp(CorrectOrdering,:)];
            EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}; EnvironmentLabelsOutputTempTemp(CorrectOrdering,:)];
            LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; LegsContractedTempTemp(CorrectOrdering,:)];
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Now add in the dropped triads
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~isempty(DroppedTriads)
            for kk = size(DroppedTriads,1):-1:1
                Loc = find(any(DroppedTriads(kk,3) == ActionTensorsTemp{oo,ss}(:,[1,2]),2),1,'first');
                if ~isempty(Loc)
                    
                    ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}(1:(Loc-1),:); DroppedTriads(kk,:);ActionTensorsTemp{oo,ss}(Loc:end,:)];
                    EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss}(1:(Loc-1),:); ForAction_EnvironmentLabelsForDroppedTriads(kk,:);EnvironmentLabelsOutputTemp{oo,ss}(Loc:end,:)];
                    LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}(1:(Loc-1),:); ForAction_DroppedTriadLegContractions(kk,:);LegsContractedTemp{oo,ss}(Loc:end,:)];
                    
                elseif ForAction_DroppedTriadNeededForBelow(kk)
                    
                    ActionTensorsTemp{oo,ss} = [ActionTensorsTemp{oo,ss}; DroppedTriads(kk,:)];
                    EnvironmentLabelsOutputTemp{oo,ss} = [EnvironmentLabelsOutputTemp{oo,ss};  ForAction_EnvironmentLabelsForDroppedTriads(kk,:)];
                    LegsContractedTemp{oo,ss} = [LegsContractedTemp{oo,ss}; ForAction_DroppedTriadLegContractions(kk,:)];
                    
                end
            end
        end
        
        
        ActionTypeTemp{oo,ss} = repmat(2,[size(ActionTensorsTemp{oo,ss},1),1]);
        for kk = 1:numel(ActionTypeTemp{oo,ss})
            if isempty(LegsContractedTemp{oo,ss}{kk})
                ActionTypeTemp{oo,ss}(kk) = 0;
            end
        end
        
        
    end %this is for the for loop.
    
    ActionTotalTensors = zeros([0,3]);
    EnvironmentLabelsOutputTotal = cell([0,1]);
    LegsContractedTotal = cell([0,1]);
    ActionTypeTotal = zeros([0,1]);
    
    for jjAll = Listing
        
        oo = Listing(1); %Contraction ordering
        ss = Listing(2); %Disjoint set number
        
        ActionTotalTensors = [ActionTotalTensors; ActionTensorsTemp{oo,ss}];
        EnvironmentLabelsOutputTotal = [EnvironmentLabelsOutputTotal; EnvironmentLabelsOutputTemp{oo,ss}];
        LegsContractedTotal = [LegsContractedTotal; LegsContractedTemp{oo,ss}];
        ActionTypeTotal = [ActionTypeTotal; ActionTypeTemp{oo,ss}];
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Simplify the network so that we don't do the same calculation multiple
    %times, this should only occur if we have multiple contraction orders
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Counter = 1;
    while Counter<size(ActionTotalTensors,1)
        
        
        
        ReplacementName = ActionTotalTensors(Counter,3);
        ContractedTensors = ActionTotalTensors(Counter,1:2);
        ContractedLegs = sort(LegsContractedTotal{Counter}(:));
        
        for kk = size(ActionTotalTensors,1):-1:(Counter+1)
            if all(ActionTotalTensors(kk,1:2) == ContractedTensors) || all(ActionTotalTensors(kk,[2,1]) == ContractedTensors)
                if isequal(ContractedLegs, sort(LegsContractedTotal{kk}(:)))
                    
                    OldName = ActionTotalTensors(kk,3);
                    EnvironmentLabelsOutputTotal{Counter} = [EnvironmentLabelsOutputTotal{Counter}, EnvironmentLabelsOutputTotal{kk}];
                    
                    ActionTotalTensors(kk,:) = [];
                    EnvironmentLabelsOutputTotal(kk) = [];
                    LegsContractedTotal(kk) = [];
                    ActionTypeTotal(kk,:) = [];
                    
                    ActionTotalTensors(ActionTotalTensors == OldName) = ReplacementName;
                end
            end
        end
        Counter = Counter+1;        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %return the traces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for kk = 1:TotalTensors(end)
        Loc = [];
        %if the tensor was traced, and if it is either one of the inputs or
        %is required for an environment output.
        if any(AtomicTensorTracing == kk) &&(any(any(ActionTotalTensors(:,[1,2]) == kk, 2)) || RequiredAtomicTensors(kk))
            Loc = find(any(ActionTotalTensors(:,[1,2]) == kk,2),1,'first');
            if isempty(Loc)
                Loc = size(ActionTotalTensors,1);
            end
            RequiredAtomicTensors(kk) = false;
            
            AtomicTensorTracingInput = kk;
            AtomicTensorTracingOutput = AtomicTensorOutputTracing(AtomicTensorTracing == kk);
        end
        
        if ~isempty(Loc)
            ActionTotalTensors(ActionTotalTensors == AtomicTensorTracingInput) = AtomicTensorTracingOutput;
            ActionTotalTensors = [ActionTotalTensors(1:(Loc-1),:); [AtomicTensorTracingInput, 0, AtomicTensorTracingOutput];ActionTotalTensors(Loc:end,:)];
            
            EnvironmentLabelsOutputTotal = [EnvironmentLabelsOutputTotal(1:(Loc-1),:); EnvironmentLabelsForAtomicTensors(kk);EnvironmentLabelsOutputTotal(Loc:end,:)];
            
            LegsContractedTotal = [LegsContractedTotal(1:(Loc-1),:); TracedLegsContractions(kk); LegsContractedTotal(Loc:end,:)];
            ActionTypeTotal = [ActionTypeTotal(1:(Loc-1),:);1;ActionTypeTotal(Loc:end,:)];
        end    
        
    end
    
    %set up the legs
    
    ActualLegsTotal = cell(size(ActionTotalTensors));
    for kk = 1:TotalTensors(end)
        Loc = find(kk>TotalTensors,1,'last');
        LegsTemp = Legs{Loc}(kk-TotalTensors(Loc),:);
        ActualLegsTotal(ActionTotalTensors == kk) = {LegsTemp};
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check that the legs required appear and create the new legs
    %worry about unmovable legs here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for ll = 1:size(ActionTotalTensors,1)
        
        if ActionTypeTotal(ll) == 2 %if contraction
            if isempty(ActualLegsTotal{ll,1})||isempty(ActualLegsTotal{ll,2})
                error('Assertion Error: Since this is a contraction then there must be a leg on both tensors')
            end
            Keep1 = ~any(repmat(ActualLegsTotal{ll,1}, [numel(LegsContractedTotal{ll}),1]) == repmat(LegsContractedTotal{ll}(:), [1,size(ActualLegsTotal{ll,1},2)]),1);
            Keep2 = ~any(repmat(ActualLegsTotal{ll,2}, [numel(LegsContractedTotal{ll}),1]) == repmat(LegsContractedTotal{ll}(:), [1,size(ActualLegsTotal{ll,2},2)]),1);
            ActualLegsTotal{ll,3} = [ActualLegsTotal{ll,1}(Keep1), ActualLegsTotal{ll,2}(Keep2)];
            
            %do a check
            Check1 = sum(repmat(ActualLegsTotal{ll,1}, [numel(LegsContractedTotal{ll}),1]) == repmat(LegsContractedTotal{ll}(:), [1,size(ActualLegsTotal{ll,1},2)]),2);
            Check2 = sum(repmat(ActualLegsTotal{ll,2}, [numel(LegsContractedTotal{ll}),1]) == repmat(LegsContractedTotal{ll}(:), [1,size(ActualLegsTotal{ll,2},2)]),2);
            if ~(all(Check1 == 1) && all(Check2 == 1))
                error('failed to Afirm, the contraction is wrong')
            end
            
            LocationsOfNewLegs = [false([ll,3]);[ActionTotalTensors((ll+1):end,1:2) == ActionTotalTensors(ll,3), false([size(ActionTotalTensors,1)-ll,1])]  ];
            
            ActualLegsTotal(LocationsOfNewLegs(:)) = repmat(ActualLegsTotal(ll,3), [sum(LocationsOfNewLegs(:)),1]);
            
        else %if trace
            if isempty(ActualLegsTotal{ll,1})
                error('Assertion Error: Since this is a trace then there must be at least two legs on the tensor')
            end
            Keep1 = ~any(repmat(ActualLegsTotal{ll,1}, [numel(LegsContractedTotal{ll}),1]) == repmat(LegsContractedTotal{ll}(:), [1,size(ActualLegsTotal{ll,1},2)]),1);
            ActualLegsTotal{ll,3} = ActualLegsTotal{ll,1}(Keep1);
            
            %do a check
            Check1 = sum(repmat(ActualLegsTotal{ll,1}, [numel(LegsContractedTotal{ll}),1]) == repmat(LegsContractedTotal{ll}(:), [1,size(ActualLegsTotal{ll,1},2)]),2);
            if ~(all(Check1 == 2))
                error('failed to Afirm, the contraction is wrong')
            end
            
            LocationsOfNewLegs = [false([ll,3]);[ActionTotalTensors((ll+1):end,1:2) == ActionTotalTensors(ll,3), false([size(ActionTotalTensors,1)-ll,1])]  ];
            
            ActualLegsTotal(LocationsOfNewLegs(:)) = repmat(ActualLegsTotal(ll,3), [sum(LocationsOfNewLegs(:)),1]);
                        
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %from this point on it doesn't matter about the specific label that we
    %are using in ActionTotalTensors therefore we will minimise the number
    %used with the exception of the atomic tensors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ActionTotalTensors = ActionTotalTensors';
    
    [~,Index,Replacements] = unique(ActionTotalTensors(ActionTotalTensors>TotalTensors(end)|ActionTotalTensors<0));
    [~,~,Index] = unique(Index);
    Replacements = Index(Replacements);
    
    ActionTotalTensors(ActionTotalTensors>TotalTensors(end)|ActionTotalTensors<0) = Replacements+TotalTensors(end);
    ActionTotalTensors = ActionTotalTensors';
    
    
    %store which locations will be used for the environment, this requires
    %both the location and the disjoint set this is from
    
    EnvironmentAccess = cell([size(LegOrder,1),1]);
    EnvironmentAccessNumbers = zeros([size(LegOrder,1)+1,1]);
    EnvironmentAccessDetails = zeros([4,0]);
    
    for nn = 1:size(EnvironmentAccess,1)
        
        EnvironmentAccessDetailsTemp = zeros([4,0]);
        
        for kk = 1:numel(EnvironmentLabelsForAtomicTensors)
            if ~isempty(EnvironmentLabelsForAtomicTensors{kk}) && ~any(AtomicTensorTracing == kk)
                Temp = EnvironmentLabelsForAtomicTensors{kk}(1,:) ==  nn;
                if any(Temp)
                    EnvironmentAccessDetailsTemp = [EnvironmentAccessDetailsTemp, [repmat([kk;0;kk], [1,sum(Temp)]);EnvironmentLabelsForAtomicTensors{kk}(2,Temp)]];
                end
            end
        end
        
        
        for kk = 1:numel(EnvironmentLabelsOutputTotal)
            if ~isempty(EnvironmentLabelsOutputTotal{kk})
                Temp = EnvironmentLabelsOutputTotal{kk}(1,:) ==  nn;
                if any(Temp)
                    EnvironmentAccessDetailsTemp = [EnvironmentAccessDetailsTemp, [repmat([ActionTotalTensors(kk,3);kk;ActionTotalTensors(kk,3)], [1,sum(Temp)]);EnvironmentLabelsOutputTotal{kk}(2,Temp)]];
                end
            end
        end
        %now I have a variable where it is a cell with a size equal to the 
        %number of outputs that we wanted, for each cell there is a double 
        %column entry where the first column is the number of the
        %contraction which generates one of the disjoint sets of this
        %system. The second column is row in ActionTotalTensors where this
        %appears. The third column repeats the label in ActionTensorsTotal 
        %that the output is set to be. The fourth column is the disjoint 
        %set that this part contributes. We will order with respect to the 
        %disjoint sets.
        
        [~,Index] = sort(EnvironmentAccessDetailsTemp(4,:),'ascend');
        EnvironmentAccessDetailsTemp = EnvironmentAccessDetailsTemp(:,Index);
        
        EnvironmentAccessDetails = [EnvironmentAccessDetails, EnvironmentAccessDetailsTemp];
        EnvironmentAccessNumbers(nn+1) = size(EnvironmentAccessDetails,2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now store if the tensor is going to be stored or removed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    KeepUnique = true([size(ActionTotalTensors,1),2]);
    
    ActionTensorsLocations = ActionTotalTensors;
    
    %now work out locations, first by changing all multiple copies of a
    %tensor into the single location
    for nn = 1:(length(TotalTensors)-1)
        RelevantLocations = (ActionTotalTensors<=TotalTensors(nn+1))&(ActionTotalTensors>TotalTensors(nn));
        ActionTensorsLocations(RelevantLocations) = nn;
        RelevantLocations = (EnvironmentAccessDetails(1,:)<=TotalTensors(nn+1))&(EnvironmentAccessDetails(2,:)>TotalTensors(nn));
        EnvironmentAccessDetails(1,RelevantLocations) = nn;
    end
    
    ActionTensorsLocations = ActionTensorsLocations-(TotalTensors(end)-(length(TotalTensors)-1))*(ActionTensorsLocations>TotalTensors(end));
    EnvironmentAccessDetails(1,:) = EnvironmentAccessDetails(1,:)-(TotalTensors(end)-(length(TotalTensors)-1))*(EnvironmentAccessDetails(1,:)>TotalTensors(end));
    
    %we need to use ActionTensorsLocations here because this way we have
    %already compressed the input tensors back into being the same as would
    %be needed when we need to work out if we still need to allocate memory
    %for them.
    KeptForKronicasAtEnd = true([max(max(ActionTensorsLocations)),1]);
    for nn = 1:max(max(ActionTensorsLocations))
        %first work out when it is output
        
        
        if nn>=numel(TotalTensors)
            OutputLoc = find(ActionTensorsLocations(:,3)==nn);
            if numel(OutputLoc) ~=1
                error('Assertion Failed: OutputLoc should have exactly 1 entry');
            end
            if sum(sum(ActionTensorsLocations(1:OutputLoc,1:2) == nn))>0
                error('Assertion Failed: OutputLoc should have been created before it was ever used as an input');
            end
            
            WorkingEnvironments = EnvironmentLabelsOutputTotal{OutputLoc};
        else
            WorkingEnvironments = EnvironmentLabelsForAtomicTensors{nn};
        end
        
        if isempty(WorkingEnvironments)
            LastColumn1 = find(ActionTensorsLocations(:,1) == nn,1, 'last');
            LastColumn2 = find(ActionTensorsLocations(:,2) == nn,1, 'last');
            
            if isempty(LastColumn1) && isempty(LastColumn2) && nn>TotalTensors(end)
                error('Assertion Failed: There should have been at least one use of this tensor since it was created (or it should be an output');
            end
            
            if isempty(LastColumn1)
                KeepUnique(LastColumn2,2) = false;
            elseif isempty(LastColumn2)
                KeepUnique(LastColumn1,1) = false;
            else
                if LastColumn1 == LastColumn2
                    error('Assertion Failed: the same tensor can''t be used twice in a single triad');
                end
                
                LastRow = max(LastColumn1, LastColumn2);
                if LastRow == LastColumn1
                    KeepUnique(LastColumn1,1) = false;
                else
                    KeepUnique(LastColumn2,2) = false;
                end
            end
            KeptForKronicasAtEnd(nn) = false;
        end
    end
    
    
    %Fix up locations to minimise the number of cells I use:
    for kk = 1:size(KeepUnique,1)
        if any(~KeepUnique(kk,:))
            UpdateOrder = sort(unique(ActionTensorsLocations(kk,~KeepUnique(kk,:))),'descend');
            
            for UpdateLoc = UpdateOrder
                ActionTensorsLocations((kk+1):end,:) = ActionTensorsLocations((kk+1):end,:) - (ActionTensorsLocations((kk+1):end,:)>UpdateLoc);
                ActionTensorsLocations(kk,3) = ActionTensorsLocations(kk,3) - (ActionTensorsLocations(kk,3)>UpdateLoc);
                EnvironmentAccessDetails(1,:) = EnvironmentAccessDetails(1,:) - (EnvironmentAccessDetails(1,:)>UpdateLoc);%&(EnvironmentAccessDetails(2,:)>kk));
            end
            
        end
    end
    
    for nn = 1:numel(EnvironmentAccess)
        EnvironmentAccess{nn} = EnvironmentAccessDetails(:,(EnvironmentAccessNumbers(nn)+1):EnvironmentAccessNumbers(nn+1));
    end
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This part is to take into account Kronica products that happen at the
    %end of the contraction process (later this may be replaced by treating
    %the tensors as a tensor product) and to construct OldTrue.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ActualLegsKronica = cell([0,3]);
    ActionTensorsKronica = zeros([0,3]);
    ActionTensorsLocationsKronica = zeros([0,3]);
    KeepUniqueKronica = zeros([0,2]);
    
    CurrentAccess = max(max(ActionTotalTensors));
    OldTrue = zeros([1,numel(EnvironmentAccess)]);
    MaxLocation = sum(KeptForKronicasAtEnd);
    
    
    OldTrueLocations = zeros([1,numel(EnvironmentAccess)]);
    EnvironmentLegs = cell([numel(EnvironmentAccess),1]);
    
    for ee = 1:numel(EnvironmentAccess)
        
        OldTrue(ee) = EnvironmentAccess{ee}(3,1);
        if OldTrue(ee)>TotalTensors(end);
            Loc = find(ActionTotalTensors(:,3) == OldTrue(ee),1,'first');
            Legs1 = ActualLegsTotal{Loc,3};
        else
            Loc = sum(OldTrue(ee)>TotalTensors);
            Legs1 = Legs{Loc}(OldTrue(ee)-TotalTensors(Loc),:);
        end
        
        OldTrueLocations(ee) = EnvironmentAccess{ee}(1,1);
        KeepUnique1 = any(EnvironmentAccessDetails(3,(EnvironmentAccessNumbers(ee)+2):end) == OldTrue(ee));
        
        for nn = 2:size(EnvironmentAccess{ee},2)
            
            %update new values
            NewLocation = EnvironmentAccess{ee}(1,nn);
            NewNumber = EnvironmentAccess{ee}(3,nn);
            KeepUnique2 = any(EnvironmentAccessDetails(3,(EnvironmentAccessNumbers(ee)+nn+1):end) == NewNumber);
            MaxLocation = MaxLocation - ~KeepUnique1 - ~KeepUnique2;
            
            %add the values to each of the Kronica parameters
            ActionTensorsKronica = [ActionTensorsKronica; [OldTrue(ee),NewNumber, CurrentAccess+1]];
            ActionTensorsLocationsKronica = [ActionTensorsLocationsKronica;[OldTrueLocations(ee), NewLocation, MaxLocation]];
            KeepUniqueKronica = [KeepUniqueKronica; [KeepUnique1, KeepUnique2]];
            
            %update states for the next iteration or for output
            CurrentAccess = CurrentAccess+1;
            OldTrue(ee) = CurrentAccess;
            OldTrueLocations(ee) = MaxLocation;
            KeepUnique1 = false;
            
            %work out the Legs
            if NewNumber>TotalTensors(end);
                Loc = find(ActionTotalTensors(:,3) == NewNumber,1,'first');
                Legs2 = ActualLegsTotal{Loc,3};
            else
                Loc = sum(NewNumber>TotalTensors);
                Legs2 = Legs{Loc}(NewNumber-TotalTensors(Loc),:);
            end
            
            %extend the Kronica
            ActualLegsKronica = [ActualLegsKronica; {Legs1, Legs2, [Legs1,Legs2]}];
            %updateLegs1 for the next iteration or for the output
            Legs1 = [Legs1,Legs2];
        end
        EnvironmentLegs{ee} = Legs1;
    end
    
    ActionTotalTensors = [ActionTotalTensors;ActionTensorsKronica];
    ActionTensorsLocations = [ActionTensorsLocations;ActionTensorsLocationsKronica];
    ActualLegsTotal = [ActualLegsTotal;ActualLegsKronica];
    DropTensor = ~[KeepUnique; KeepUniqueKronica];
    
    ActionTypeTotal = [ActionTypeTotal;zeros([size(ActionTensorsKronica,1),1])];
    LegsContractedTotal = [LegsContractedTotal;cell([size(ActionTensorsKronica,1),1])];
end


function [OldTrueChargeSidesInt, OldTrueBraiding, OldTrueBraidingDirection, OldTrueStructure, OldTrueChargeDirections,...
        ActionTensorsPermute,ActionNumberContractedLegs,OldTrueChargeSides,OldTruePermute,ActionTensorsStructure,ActionTensorsChargeDirections] = ...
    Broad_5_ComputingOldTrue(LegOrder, EnvironmentLegs, StructureTensorDetails, ChargeDirectionsTensorDetails, ActualLegsTotal, LegsContractedTotal, ...
    ActionTotalTensors, ActionTypeTotal, ActionTensorsLocations, OldTrue, OldTrueLocations, TotalTensors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%STEP 5 - BROAD: Compute the operations we have to do on the output tensors
%assuming that there are no symmetries. We also do the alterations due to
%permutations here as well
%
%Inputs - LegOrder [1]
%       - EnvironmentLegs [4]
%       - StructureTensorDetails [Input]
%       - ChargeDirectionsTensorDetails [Input]
%       - ActualLegsTotal [4]
%       - LegsContractedTotal [4]
%       - ActionTotalTensors [4]
%       - ActionTypeTotal [4]
%
%Outputs 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %here I set the parameters for the output tensors (assuming they are
    %not symmetric, if they are symmetric then I will update this later).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OldTrueChargeSides = cell([size(LegOrder,1),1]);
    PermuteFinal = cell([size(LegOrder,1),1]);
    OldTrueChargeSidesInt = cell([size(LegOrder,1),1]);
    OldTrueBraiding = cell([size(LegOrder,1),1]);
    OldTrueBraidingDirection = cell([size(LegOrder,1),1]);
    
    for kk = 1:size(LegOrder,1)
        
        A = EnvironmentLegs{kk};
        B = LegOrder{kk,2};
        
        %compute the permutation required
        [~,IndexA] = sort(A,'ascend');
        [~,IndexB] = sort(B,'ascend');
        [~,IndexB2] = sort(IndexB,'ascend'); %note that this is the inverse ordering of leg ordering B
        
        PermuteFinal{kk}  = IndexA(IndexB2);
        [Check,FinalLocationAccess] = sort(OldTrueLocations,'ascend');
        
        if ~isequal(Check, 1:length(Check))
            error('Affirmation Failed: I expected only the final tensors to be left by now');
        end
        OldTrueChargeSides{kk} = LegOrder{kk,4};
        OldTrueChargeSidesInt{kk} = LegOrder{kk,5};
        OldTrueBraiding{kk} = [];
        OldTrueBraidingDirection{kk} = [];
    end
    
    OldTruePermute = PermuteFinal;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now to work out how I should be permuting for ActionTensorsTotal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ActionTensorsPermute = cell(size(ActionTensorsLocations));
    ActionNumberContractedLegs = zeros([size(ActionTensorsLocations,1),1]);
    
    
    for kk = 1:size(ActualLegsTotal,1)
        ActionNumberContractedLegs(kk) = numel(LegsContractedTotal{kk});
        
        if ActionTypeTotal(kk) ==2 %if a contraction
            A = any(repmat(LegsContractedTotal{kk}(:),[1,size(ActualLegsTotal{kk,1},2)]) == repmat(ActualLegsTotal{kk,1},[numel(LegsContractedTotal{kk}),1]), 1);
            [AOrderedLogic,IndexA] = sort(A,'ascend');
            [CheckA, TrueIndexA] = sort(ActualLegsTotal{kk,1}(A),'ascend');
            TempIndexA = IndexA(AOrderedLogic);
            IndexA(AOrderedLogic) = TempIndexA(TrueIndexA);
            
            B = any(repmat(LegsContractedTotal{kk}(:),[1,size(ActualLegsTotal{kk,2},2)]) == repmat(ActualLegsTotal{kk,2},[numel(LegsContractedTotal{kk}),1]), 1);
            [BOrderedLogic,IndexB] = sort(B,'descend');
            [CheckB, TrueIndexB] = sort(ActualLegsTotal{kk,2}(B),'ascend');
            TempIndexB = IndexB(BOrderedLogic);
            IndexB(BOrderedLogic) = TempIndexB(TrueIndexB);
            
            if ~isequal(CheckA, CheckB)
                error('Affirmation Failed: the contraction permutation turned out wrong')
            end
            
            Permute1 = IndexA;
            Permute2 = IndexB;
            
            Permute3 = [];
            
            
        elseif ActionTypeTotal(kk) ==1 %if a trace
            A = any(repmat(LegsContractedTotal{kk}(:),[1,size(ActualLegsTotal{kk,1},2)]) == repmat(ActualLegsTotal{kk,1},[numel(LegsContractedTotal{kk}),1]), 1);
            [AOrderedLogic,IndexA] = sort(A,'descend');
            [Check, TrueIndexA] = sort(ActualLegsTotal{kk,1}(A),'ascend');
            TempIndexA = IndexA(AOrderedLogic);
            if ~isequal(Check(1:2:end), Check(2:2:end))
                error('Affirmation Failed: the trace permutation turned out wrong')
            end
            IndexA(AOrderedLogic) = TempIndexA(TrueIndexA([1:2:end,2:2:end]));
            
            Permute1 = IndexA;
            
            Permute2 = [];
            Permute3 = [];
            
        else %ActionTypesTotal(kk) ==0 %if a kronica
            Permute1 = [];
            Permute2 = [];
            Permute3 = [];
        end
        
        ActionTensorsPermute(kk,1:3) = {Permute1,Permute2,Permute3};
    end
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now work out Structure and ChargeDirections for all resultant
    %states, this is simple because we are assuming it is in standard form.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ActionTensorsStructure = cell(size(ActionTotalTensors));
    ActionTensorsChargeDirections = cell(size(ActionTotalTensors));
    
    for nn = 1:(numel(TotalTensors)-1)
        Temp = ActionTotalTensors>TotalTensors(nn)&ActionTotalTensors<=TotalTensors(nn+1);
        ActionTensorsStructure(Temp) = repmat(StructureTensorDetails(nn), [sum(sum(Temp)),1]);
        ActionTensorsChargeDirections(Temp) = repmat(ChargeDirectionsTensorDetails(nn), [sum(sum(Temp)),1]);
    end
    
    for nn = 1:size(ActionTotalTensors,1)
        
        %Construct Structure
        NumberSizes = size(ActualLegsTotal{nn,3},2);
        StructureUsed = SymTensor.GenerateStructure([NumberSizes,1]);
        
        ActionTensorsStructure(ActionTotalTensors(:)==ActionTotalTensors(nn,3)) = repmat({StructureUsed},[sum(ActionTotalTensors(:)==ActionTotalTensors(nn,3)),1]);
        
        ChargeDirectionsUsed = ones(size(ActionTensorsStructure{nn,3}));
        %if there are more then one leg on the nn th output then work out
        %the associated charge direction of this output leg from inputs
        if NumberSizes > 1
            Charges = zeros(size(ActualLegsTotal{nn,3}));
            for ll = 1:NumberSizes
                LegName = ActualLegsTotal{nn,3}(ll);
                OldNumber1 = find(ActualLegsTotal{nn,1}==LegName);
                if ActionTypeTotal(nn)==1
                    OldNumber2 = [];
                    if numel(OldNumber1) ~=1
                        error('Affirmation Failed: We should have exactly 1 leg lining up in this manner')
                    end
                else
                    OldNumber2 = find(ActualLegsTotal{nn,2}==LegName);
                    if numel(OldNumber1)+numel(OldNumber2) ~=1
                        error('Affirmation Failed: We should have exactly 1 leg lining up in this manner')
                    end
                end
                
                if numel(OldNumber1) == 1
                    LegNumber = OldNumber1;
                    Charges(ll) = ActionTensorsChargeDirections{nn,1}(ActionTensorsStructure{nn,1}==-LegNumber);
                end
                
                if numel(OldNumber2) == 1;
                    LegNumber = OldNumber2;
                    Charges(ll) = ActionTensorsChargeDirections{nn,2}(ActionTensorsStructure{nn,2}==-LegNumber);
                end
                ChargeDirectionsUsed(ActionTensorsStructure{nn,3} == -ll) = Charges(ll);
            end
            
        end
        %else we realise that the vector has to be a trivial irrep and
        %therefore will be self duel so the charge direction doesn't
        %matter.
        
        ActionTensorsChargeDirections(ActionTotalTensors(:)==ActionTotalTensors(nn,3)) = repmat({ChargeDirectionsUsed},[sum(ActionTotalTensors(:)==ActionTotalTensors(nn,3)),1]);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now work out Structure and ChargeDirections for all output tensors
    %(the OldTrue tensors)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    OldTrueStructure = cell([1,size(OldTrue,2)]);
    OldTrueChargeDirections = cell([1,size(OldTrue,2)]);
    OldTrueStoredLocations = cell([1,size(OldTrue,2)]);
    
    for nn = 1:size(OldTrue,2)
        Temp = ActionTotalTensors(:,3) == OldTrue(nn);
        if sum(sum(Temp))~=1
            error('Affirmation Failed: There should only be one copy of this produced')
        end
        StructureUsed = LegOrder{nn,3};
        ChargeDirectionsUsed = ones(size(StructureUsed));
        
        %if there is more then one leg then take the structure from the
        %precomputed ones
        if ~isempty(StructureUsed)
            
            EffectiveChargeDirections = ActionTensorsChargeDirections{ActionTotalTensors(:,3)==OldTrue(nn),3};
            EffectiveStructure = ActionTensorsStructure{ActionTotalTensors(:,3)==OldTrue(nn),3};
            
            %work out charge directions for output legs
            for ll = 1:(size(StructureUsed,2)+1)
                ChargeDirectionsUsed(StructureUsed==-ll) = EffectiveChargeDirections(EffectiveStructure == -OldTruePermute{nn}(ll));
            end
            
            %work out charge directions for internal legs
            for ll = 1:(size(StructureUsed,2)-1)
                ChargeDirectionsUsed(StructureUsed==ll) = OldTrueChargeSidesInt{nn}(ll);
            end
        end
        OldTrueStructure{nn} = StructureUsed;
        OldTrueChargeDirections{nn} = ChargeDirectionsUsed;
        OldTrueStoredLocations{nn} = SymTensor.GenerateStoredLocations(StructureUsed);
    end
end