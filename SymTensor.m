classdef SymTensor
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess = 'public', SetAccess = 'public')
        TensorEntries %The parts of the tensors 
        TensorEntriesSizes %the sizes of the tensors
        
        
        Structure %How to combine the tripoints and endpoints, note that we start counting from zero.
        ChargeDirections %what directions the tripoints are going in, +1 = Out of Core, -1 = Into core
        
        ChargeLabelsExternal %the external charges, they change only due to braidings
        ChargeLabelsInternal %the internal charges, they change with internal restructuring
        MultiplicitiesInternal %the internal multiplicities, they can change with internal restructuring
        
        Braidings %a list of the charges braided, this is always the leftmost of CyclicLegs.
        %note that this doesn't reference the output leg as a label but is an internal thing.
        BraidingDirections %directions for the braidings, +1 = line going upright is above, -1 = line going upleft is above, 
        %note that this is regardless of where it is which means it can
        %have fiferent meanings and actions on top side and bottom
        StoredLocations %quick reference for the ordering at the output.
        CapType %tells us what the cap and cups are like +1 means going to right, -1 means going to the left
        
        
        ChargeLegDimensions % a list of dimensions associated to each charge on each leg. indexes are charges, then legs
        ChargeSide %a list of which side of the block the leg is, +1 = top, -1 = bottom
        ChargeSideInt %a list of which side of the vertex it directs into the leg is, -1 = top, +1 = bottom; 
                      %NOTE: this is different to ChargeSide only because
                      %we want these to be +1 if on top half (going into
                      %bottom of vertex) and -1 if on bottom half (going
                      %into top of vertex).
        SymHandle %the Symmetry associated
        
        FlagIsNotSymmetric %If this tensors has only used the trival irrep
        FlagIsANumber %If this tensors is actually a number
        FlagIsAShape = false; %if this tensor is actually a shape
        
        %FreePermuteTracker
        
        CyclicLegs %this is the leg ordering after the braiding
        CyclicLegsBraided %this is the leg ordering when we change to have top then bottom, note that this is when we have braidings in addition to this.
        
    end
    
    
    properties(GetAccess = 'public', SetAccess = 'private')
        
    end
    
    methods(Access = 'public', Static = true)
        
        function OutTensor = CreateTensor(Sym, TypeMatrix, ChargeDimensionsBottom, ChargeDimensionsTop, NumberOfLegsBottom, NumberOfLegsTop)
            
            if nargin<4
                ChargeDimensionsTop = ChargeDimensionsBottom;
            end
            if nargin<5
                NumberOfLegsBottom = size(ChargeDimensionsBottom,2);
            end
            
            if nargin<6
                NumberOfLegsTop = size(ChargeDimensionsTop,2);
            end
            
            
            %now check that the number of dimensions is consistant
            
            TempDim = Sym.getDimExternal();
            
            if size(ChargeDimensionsBottom,1) ~= size(ChargeDimensionsTop,1)
                error('error: ChargeDimensionsBottom and ChargeDimensionsTop should agree')
            end
            
            if length(size(ChargeDimensionsBottom))>3
                error('error: ChargeDimensionsBottom should be a vector, not 3D')
            end
            if length(size(ChargeDimensionsTop))>3
                error('error: ChargeDimensionsTop should be a vector, not 3D')
            end
            
            
            
            if ~(size(ChargeDimensionsBottom,2)==1)
                %Then either we are assigning seperately
                if size(ChargeDimensionsBottom,2)==sum(abs(NumberOfLegsBottom))
                    
                    if (size(ChargeDimensionsBottom,1)~=size(TempDim,1))
                        error('error: ChargeDimensionsBottom should have the same number of entries as the number of Charges in the symmetry being used')
                    end
                    
                    FlagAssignSeperateLegsBottom = true;
                    
                elseif size(ChargeDimensionsBottom,1)==size(TempDim,1) && size(ChargeDimensionsBottom,2) == 1
                    FlagAssignSeperateLegsBottom = false;
                else
                    error('error: ChargeDimensionsBottom should be a vector, not a matrix')
                end
            else
                if (size(ChargeDimensionsBottom,1)~=size(TempDim,1))
                    error('error: ChargeDimensions should have the same number of entries as the number of Charges in the symmetry being used')
                end
                FlagAssignSeperateLegsBottom = false;
            end
            
            
            
            if any(ChargeDimensionsBottom~=round(ChargeDimensionsBottom))
                error('error: the degeneracy dimension must be an integer (specifically positive non-negative)')
            end
            
            if any(ChargeDimensionsBottom~=real(ChargeDimensionsBottom))
                error('error: the degeneray dimension must be a real integer (specifically non-negative)')
            end
            
            if any(ChargeDimensionsBottom < 0)
                error('error: The degeneracy dimension must be zero or a positive number')
            end
            
            
            if ~(size(ChargeDimensionsTop,2)==1)
                %Then either we are assigning seperately
                if size(ChargeDimensionsTop,2)==sum(abs(NumberOfLegsTop))
                    
                    if (size(ChargeDimensionsTop,1)~=size(TempDim,1))
                        error('error: ChargeDimensionsTop should have the same number of entries as the number of Charges in the symmetry being used')
                    end
                    
                    FlagAssignSeperateLegsTop = true;
                    
                elseif size(ChargeDimensionsTop,1)==size(TempDim,1) && size(ChargeDimensionsTop,2) == 1
                    FlagAssignSeperateLegsTop = false;
                else
                    error('error: ChargeDimensionsTop should be a vector, not a matrix')
                end
            else
                if (size(ChargeDimensionsTop,1)~=size(TempDim,1))
                    error('error: ChargeDimensionsTop should have the same number of entries as the number of Charges in the symmetry being used')
                end
                FlagAssignSeperateLegsTop = false;
            end
            
            
            
            if any(ChargeDimensionsTop~=round(ChargeDimensionsTop))
                error('error: the degeneracy dimension must be an integer (specifically positive non-negative)')
            end
            
            if any(ChargeDimensionsTop~=real(ChargeDimensionsTop))
                error('error: the degeneray dimension must be a real integer (specifically non-negative)')
            end
            
            if any(ChargeDimensionsTop < 0)
                error('error: The degeneracy dimension must be zero or a positive number')
            end
            
            
            
            
            %now begin
            
            %first work out if there are multiplicities, and if it is
            %abelian, self dualness is in the next part when we work out
            %the legs properly, Braiding is unimportant
            
            FlagNonAbelian = IsNonAbelian(Sym);
            FlagMultiplicities = IsMultiplicities(Sym);
            
            
            Fusion2 = getFusion2(Sym);
                
            
            %now create list of Allowed charges
            
            AllowedSingleCharges = TempDim(:,1)';
            AllowedSingleCharges = AllowedSingleCharges(any([ChargeDimensionsBottom,ChargeDimensionsTop]~=0,2));
            
            InverseIrrep = Sym.getInverseIrrep;
            
            %now work out the legs;           
            
            TotalLegsBottom = sum(abs(NumberOfLegsBottom));
            TotalLegsTop = sum(abs(NumberOfLegsTop));
            TotalLegs = sum(abs(NumberOfLegsBottom))+sum(abs(NumberOfLegsTop));
            ChargeDirectionsUsedBottom = ones([size(TempDim,1),TotalLegsBottom]);
            ChargeDirectionsUsedTop = ones([size(TempDim,1),TotalLegsTop]);
            
            %if ~IsSelfDual(Sym)
                ChargeDirectionsUsedBottom = [];
                for jj = 1:size(NumberOfLegsBottom,1)
                    ChargeDirectionsUsedBottom = [ChargeDirectionsUsedBottom,-sign(NumberOfLegsBottom(jj))*ones([1,abs(NumberOfLegsBottom(jj))])];
                end
                
            %end
            
            
            %if ~IsSelfDual(Sym)
                ChargeDirectionsUsedTop = [];
                for jj = 1:size(NumberOfLegsTop,1)
                    ChargeDirectionsUsedTop = [ChargeDirectionsUsedTop,sign(NumberOfLegsTop(jj))*ones([1,abs(NumberOfLegsTop(jj))])];
                end
                
            %end
            
            
            
            %now work out the combinations allowed
            

            
            if TotalLegsBottom>1
                
                AllowedSingleCharges = TempDim(:,1)';
                AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsBottom(:,1)~=0)';
                
                ChargeLabelsExternalBottom = AllowedSingleCharges;
            
                if FlagAssignSeperateLegsBottom
                
                    for jj = 2:TotalLegsBottom
                        
                        AllowedSingleCharges = TempDim(:,1)';
                        AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsBottom(:,jj)~=0)';
                        
                        
                        ChargeLabelsExternalBottom = [repmat(ChargeLabelsExternalBottom, [size(AllowedSingleCharges,1), 1]),...
                            reshape(repmat(reshape(AllowedSingleCharges,[1,numel(AllowedSingleCharges)]), [size(ChargeLabelsExternalBottom,1),1]),[size(ChargeLabelsExternalBottom,1)*size(AllowedSingleCharges,1),1])];
                        
                    end
                    
                else
                    
                    for jj = 2:TotalLegsBottom
                        
                        AllowedSingleCharges = TempDim(:,1)';
                        AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsBottom(:,1)~=0)';
                        
                        
                        ChargeLabelsExternalBottom = [repmat(ChargeLabelsExternalBottom, [size(AllowedSingleCharges,1), 1]),...
                            reshape(repmat(reshape(AllowedSingleCharges,[1,numel(AllowedSingleCharges)]), [size(ChargeLabelsExternalBottom,1),1]),[size(ChargeLabelsExternalBottom,1)*size(AllowedSingleCharges,1),1])];
                        
                    end
                end
            
            ChargeLabelsInternalBottom = zeros([size(ChargeLabelsExternalBottom,1),size(ChargeLabelsExternalBottom,2)-1]);
            MultiplicitiesInternalBottom = ones([size(ChargeLabelsExternalBottom,1),size(ChargeLabelsExternalBottom,2)-1]);
            
                
                %in the abelian case we know that the sum is fixed as one
                %output, we want to store that as we are using the Symmetry
                %to fix the dimension, may add an additional option later
                %where we store if we have an abelian model even if it
                %isn't know, however this will be important in cases where
                %we have braiding where knowing the internal structure is a
                %must, so it would have to be mutually exclusive to the
                %symmetry containing braiding in general.
                
                
            %first one is special:
            %FIXHERE need to take case of one leg each side into
            %account
            
            if ChargeDirectionsUsedBottom(1) == +1
                Charge1 = InverseIrrep(ChargeLabelsExternalBottom(:,1));
            else
                Charge1 = ChargeLabelsExternalBottom(:,1);
            end
            
            
            if ChargeDirectionsUsedBottom(2) == +1
                Charge2 = InverseIrrep(ChargeLabelsExternalBottom(:,2));
            else
                Charge2 = ChargeLabelsExternalBottom(:,2);
            end
             
            ChargeLabelsInternalBottom = zeros([size(ChargeLabelsExternalBottom,1), size(ChargeLabelsExternalBottom,2)-1]);
            MultiplicitiesInternalBottom = zeros([size(ChargeLabelsExternalBottom,1), size(ChargeLabelsExternalBottom,2)-1]);
             
            %first one
            [TempChargeLabelsInternal, MultList, Multiplicities] = Sym.FuseChargeList(Charge1,Charge2);
            
            ChargeLabelsExternalBottom = ChargeLabelsExternalBottom(MultList,:);
            ChargeLabelsInternalBottom = ChargeLabelsInternalBottom(MultList,:);
            MultiplicitiesInternalBottom = MultiplicitiesInternalBottom(MultList,:);
            
            ChargeLabelsInternalBottom(:,1) = TempChargeLabelsInternal;
            MultiplicitiesInternalBottom(:,1) = Multiplicities;
            
            for kk = 2:size(ChargeLabelsInternalBottom,2)
                if ChargeDirectionsUsedBottom(kk+1) == +1
                    [TempChargeLabelsInternal, MultList, Multiplicities] = Sym.FuseChargeList(ChargeLabelsInternalBottom(:,kk-1),InverseIrrep(ChargeLabelsExternalBottom(:,kk+1)));
                else
                    [TempChargeLabelsInternal, MultList, Multiplicities] = Sym.FuseChargeList(ChargeLabelsInternalBottom(:,kk-1),ChargeLabelsExternalBottom(:,kk+1));
                end
                
                ChargeLabelsExternalBottom = ChargeLabelsExternalBottom(MultList,:);
                ChargeLabelsInternalBottom = ChargeLabelsInternalBottom(MultList,:);
                MultiplicitiesInternalBottom = MultiplicitiesInternalBottom(MultList,:);
                
                ChargeLabelsInternalBottom(:,kk) = TempChargeLabelsInternal;
                MultiplicitiesInternalBottom(:,kk) = Multiplicities;
                
            end
            
            
            elseif TotalLegsBottom == 1
                
                AllowedSingleCharges = TempDim(:,1)';
                AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsBottom(:,1)~=0);
                
                
                if ChargeDirectionsUsedBottom(1) == +1
                    ChargeLabelsExternalBottom = InverseIrrep(AllowedSingleCharges');
                else
                    ChargeLabelsExternalBottom = AllowedSingleCharges';
                end
                ChargeLabelsInternalBottom = ones([size(ChargeLabelsExternalBottom,1),0]);
                MultiplicitiesInternalBottom = ones([size(ChargeLabelsExternalBottom,1),0]);
            
            else %if TotalLegsBottom == 0
                
                ChargeLabelsExternalBottom = ones([1,0]);
                ChargeLabelsInternalBottom = ones([1,0]);
                MultiplicitiesInternalBottom = ones([1,0]);
                
            end

            
            
            %now work out the combinations allowed
            

            if TotalLegsBottom>1
                UniqueBottom = ChargeLabelsInternalBottom(:,end);
            elseif TotalLegsBottom == 1
                UniqueBottom = ChargeLabelsExternalBottom(:,end);
            else %if TotalLegsBottom == 0
                UniqueBottom = TrivialIrrep;
            end
            
            
            
            
            if TotalLegsTop>1
                
                    
                AllowedSingleCharges = TempDim(:,1)';
                AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsTop(:,1)~=0)';
                
                ChargeLabelsExternalTop = AllowedSingleCharges;
                
                if FlagAssignSeperateLegsBottom
                    
                    for jj = 2:TotalLegsTop
                        
                        AllowedSingleCharges = TempDim(:,1)';
                        AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsTop(:,jj)~=0)';
                        
                        ChargeLabelsExternalTop = [repmat(ChargeLabelsExternalTop, [size(AllowedSingleCharges,1), 1]),...
                            reshape(repmat(reshape(AllowedSingleCharges,[1,numel(AllowedSingleCharges)]), [size(ChargeLabelsExternalTop,1),1]),[size(ChargeLabelsExternalTop,1)*size(AllowedSingleCharges,1),1])];
                        
                    end
                    
                else
                    
                    for jj = 2:TotalLegsTop
                        
                        AllowedSingleCharges = TempDim(:,1)';
                        AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsTop(:,1)~=0)';
                        
                        ChargeLabelsExternalTop = [repmat(ChargeLabelsExternalTop, [size(AllowedSingleCharges,1), 1]),...
                            reshape(repmat(reshape(AllowedSingleCharges,[1,numel(AllowedSingleCharges)]), [size(ChargeLabelsExternalTop,1),1]),[size(ChargeLabelsExternalTop,1)*size(AllowedSingleCharges,1),1])];
                        
                    end
                    
                end
            ChargeLabelsInternalTop = zeros([size(ChargeLabelsExternalTop,1),size(ChargeLabelsExternalTop,2)-1]);
            MultiplicitiesInternalTop = ones([size(ChargeLabelsExternalTop,1),size(ChargeLabelsExternalTop,2)-1]);
            
            
            if ChargeDirectionsUsedTop(1) == -1
                Charge1 = InverseIrrep(ChargeLabelsExternalTop(:,1));
            else
                Charge1 = ChargeLabelsExternalTop(:,1);
            end
            
            
            if ChargeDirectionsUsedTop(2) == -1
                Charge2 = InverseIrrep(ChargeLabelsExternalTop(:,2));
            else
                Charge2 = ChargeLabelsExternalTop(:,2);
            end
             
             
            %first one
            [TempChargeLabelsInternal, MultList, Multiplicities] = Sym.FuseChargeList(Charge1,Charge2);
            
            ChargeLabelsExternalTop = ChargeLabelsExternalTop(MultList,:);
            ChargeLabelsInternalTop = ChargeLabelsInternalTop(MultList,:);
            MultiplicitiesInternalTop = MultiplicitiesInternalTop(MultList,:);
            
            ChargeLabelsInternalTop(:,1) = TempChargeLabelsInternal;
            MultiplicitiesInternalTop(:,1) = Multiplicities;
            
            for kk = 2:size(ChargeLabelsInternalTop,2)
                if ChargeDirectionsUsedTop(kk+1) == -1
                    [TempChargeLabelsInternal, MultList, Multiplicities] = Sym.FuseChargeList(ChargeLabelsInternalTop(:,kk-1),InverseIrrep(ChargeLabelsExternalTop(:,kk+1)));
                else
                    [TempChargeLabelsInternal, MultList, Multiplicities] = Sym.FuseChargeList(ChargeLabelsInternalTop(:,kk-1),ChargeLabelsExternalTop(:,kk+1));
                end
                
                ChargeLabelsExternalTop = ChargeLabelsExternalTop(MultList,:);
                ChargeLabelsInternalTop = ChargeLabelsInternalTop(MultList,:);
                MultiplicitiesInternalTop = MultiplicitiesInternalTop(MultList,:);
                
                ChargeLabelsInternalTop(:,kk) = TempChargeLabelsInternal;
                MultiplicitiesInternalTop(:,kk) = Multiplicities;
                
            end
            
            elseif TotalLegsTop == 1
                
                AllowedSingleCharges = TempDim(:,1)';
                AllowedSingleCharges = AllowedSingleCharges(ChargeDimensionsTop(:,1)~=0);
                
                if ChargeDirectionsUsedTop(1) == -1
                    ChargeLabelsExternalTop = InverseIrrep(AllowedSingleCharges');
                else
                    ChargeLabelsExternalTop = AllowedSingleCharges';
                end
                ChargeLabelsInternalTop = ones([size(ChargeLabelsExternalTop,1),0]);
                MultiplicitiesInternalTop = ones([size(ChargeLabelsExternalTop,1),0]);
            
            else 
                
                ChargeLabelsExternalTop = ones([1,0]);
                ChargeLabelsInternalTop = ones([1,0]);
                MultiplicitiesInternalTop = ones([1,0]);
                
            end

            %now work out the unique weights
            
            if TotalLegsTop>1
                UniqueTop = ChargeLabelsInternalTop(:,end);
            elseif TotalLegsTop == 1
                UniqueTop = ChargeLabelsExternalTop(:,end);
            else %if TotalLegsTop == 0
                UniqueTop = TrivialIrrep;
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %secondly we need to convert this list into a list of Matrix
            %sizes with dimensions
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            BottomEntries = [ChargeLabelsExternalBottom,ChargeLabelsInternalBottom, MultiplicitiesInternalBottom];
            TopEntries = [ChargeLabelsExternalTop,ChargeLabelsInternalTop, MultiplicitiesInternalTop];
            
            [UniqueEntries, ~, UniqueLabels] = unique([UniqueBottom;UniqueTop]);
            
            UniqueLabelsBottom = UniqueLabels(1:size(UniqueBottom,1),1);
            UniqueLabelsTop = UniqueLabels((size(UniqueBottom,1)+1):end,1);
            
            
            UniqueEntriesBottom = cell([size(UniqueEntries,1),1]);
            UniqueEntriesTop = cell([size(UniqueEntries,1),1]);
            KeepUnique = true([size(UniqueEntries,1),1]);
            
            
            ChargeDimensionsBottom = ChargeDimensionsBottom';
            ChargeDimensionsTop = ChargeDimensionsTop';
            
            for jj = 1:size(UniqueEntries,1)
                
                UniqueEntriesBottom{jj} = BottomEntries(UniqueLabelsBottom == jj,:);
                UniqueEntriesTop{jj} = TopEntries(UniqueLabelsTop == jj,:);
                
                    if TotalLegsBottom>0
                        UniqueDimensionsBottom{jj} = ChargeDimensionsBottom(1,UniqueEntriesBottom{jj}(:,1))';
                        if FlagAssignSeperateLegsBottom
                            for kk = 2:TotalLegsBottom
                                UniqueDimensionsBottom{jj} = [UniqueDimensionsBottom{jj}, ChargeDimensionsBottom(kk,UniqueEntriesBottom{jj}(:,kk))'];
                            end
                        else
                            for kk = 2:TotalLegsBottom
                                UniqueDimensionsBottom{jj} = [UniqueDimensionsBottom{jj}, ChargeDimensionsBottom(1,UniqueEntriesBottom{jj}(:,kk))'];
                            end
                        end
                    else
                        UniqueDimensionsBottom{jj} = zeros([1,0]);
                    end
                
                    if TotalLegsTop>0
                        
                        UniqueDimensionsTop{jj} = ChargeDimensionsTop(1,UniqueEntriesTop{jj}(:,1))';
                        if FlagAssignSeperateLegsTop
                            for kk = 2:TotalLegsTop
                                UniqueDimensionsTop{jj} = [UniqueDimensionsTop{jj}, ChargeDimensionsTop(kk,UniqueEntriesTop{jj}(:,kk))'];
                            end
                        else
                            for kk = 2:TotalLegsTop
                                UniqueDimensionsTop{jj} = [UniqueDimensionsTop{jj}, ChargeDimensionsTop(1,UniqueEntriesTop{jj}(:,kk))'];
                            end
                        end
                        
                    else
                        UniqueDimensionsTop{jj} = zeros([1,0]);
                    end
                    
                    KeepUnique(jj) = ~isempty(UniqueEntriesBottom{jj})&&~isempty(UniqueEntriesTop{jj});
                
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Thirdly, create the matricies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                sizes(1) = sum(prod(UniqueDimensionsBottom{jj},2));
                sizes(2) = sum(prod(UniqueDimensionsTop{jj},2));
                
                switch TypeMatrix
                    
                    case 1 %random Unitary
                        Temp = randn(max(sizes));

                        %here we make sure it is symmetric
                        Temp = 0.5*(Temp + Temp');
                        
                        %here we make sure it is traceless
                        Temp(max(sizes),max(sizes)) = -sum(diag(Temp))+Temp(max(sizes),max(sizes));
                        
                        %now change it to a unitary square matrix
                        Temp = expm(1i*Temp);
                        
                        
                        %now get the size correct (note that this is just a random projection
                        %operator occuring)
                        Matrix{jj} = Temp(1:sizes(1), 1:sizes(2));%/sqrt(Sym.Dim(jj,2));
                        
                    case 2 %identity
                        Matrix{jj} = eye(sizes);
                        
                    case 3 % ones
                        Matrix{jj} = ones(sizes);
                        

                        
                end
            end
            
            
            Matrix = Matrix(KeepUnique);
            UniqueEntriesBottom = UniqueEntriesBottom(KeepUnique);
            UniqueEntriesTop = UniqueEntriesTop(KeepUnique);
            UniqueDimensionsBottom = UniqueDimensionsBottom(KeepUnique);
            UniqueDimensionsTop = UniqueDimensionsTop(KeepUnique);
            UniqueEntries = UniqueEntries(KeepUnique);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Finally use SymTen2Mat to create the matrix
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            Legs = [TotalLegsBottom,TotalLegsTop];
            
            
            if TotalLegsBottom*TotalLegsTop >0
                Structure = zeros([2,TotalLegsBottom+TotalLegsTop-1]);
                ChargeDirections = ones([2,TotalLegsBottom+TotalLegsTop-1]);
                
                if TotalLegsBottom == 1
                    Structure(1,1) = -1;
                    ChargeDirections(1,1) = ChargeDirectionsUsedBottom(1);
                else
                    Structure(1,1) = TotalLegsBottom-1;
                    ChargeDirections(1,1) = -1;
                    Structure(2,TotalLegsBottom) = -TotalLegsBottom;
                    ChargeDirections(2,TotalLegsBottom) = ChargeDirectionsUsedBottom(TotalLegsBottom);
                end
                
                for jj = 2:(TotalLegsBottom-1)
                    Structure(1,TotalLegsBottom-jj+2) = TotalLegsBottom-jj;
                    ChargeDirections(1,TotalLegsBottom-jj+2) = -1;
                    Structure(2,TotalLegsBottom-jj+1) = -(TotalLegsBottom-jj+1);
                    ChargeDirections(2,TotalLegsBottom-jj+1) = ChargeDirectionsUsedBottom(TotalLegsBottom-jj+1);
                end
                
                if TotalLegsBottom>=2
                    Structure(1,2) = -1;
                    ChargeDirections(1,2) = ChargeDirectionsUsedBottom(1);
                end
                
                if TotalLegsTop == 1
                    Structure(2,1) = -TotalLegsBottom-1;
                    ChargeDirections(2,1) = ChargeDirectionsUsedTop(1);
                else
                    Structure(2,1) = TotalLegsBottom+TotalLegsTop-2;
                    Structure(2,TotalLegsBottom+TotalLegsTop-1) = -(TotalLegsBottom+TotalLegsTop);
                    ChargeDirections(2,TotalLegsBottom+TotalLegsTop-1) = ChargeDirectionsUsedTop(TotalLegsTop);
                end
                
                for jj = 2:(TotalLegsTop-1)
                    Structure(1,TotalLegsTop-jj+2+TotalLegsBottom-1) = TotalLegsTop-jj-1+TotalLegsBottom;
                    Structure(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = -(TotalLegsTop-jj+1+TotalLegsBottom);
                    ChargeDirections(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = ChargeDirectionsUsedTop(TotalLegsTop-jj+1);
                end
                
                if TotalLegsTop>=2
                    Structure(1,TotalLegsBottom+1) = -(TotalLegsBottom+1);
                    ChargeDirections(1,TotalLegsBottom+1) = ChargeDirectionsUsedTop(1);
                end
                
            else
                
                if TotalLegsBottom == 0
                    
                    
                    
                else %if TotalLegsTop == 0
                    
                end
            
            end
            
            Data = {UniqueEntries, UniqueDimensionsBottom, UniqueDimensionsTop, Legs, UniqueEntriesBottom, UniqueEntriesTop,...
                ChargeDirections, [],[],Sym, Structure,0};
            
            OutTensor = SymTensor.SymMat2Ten(Matrix,Data);
        end
        
        function OutTensor = CreateIdentity(Sym, ChargeDimensionsBottom, ChargeDimensionsTop, Legs)
            
            NumberOfLegsBottom = Legs(1);
            NumberOfLegsTop = Legs(2);
            
            OutTensor = SymTensor.CreateTensor(Sym,2,ChargeDimensionsBottom, ChargeDimensionsTop, NumberOfLegsBottom, NumberOfLegsTop);
            
        end
        
        function OutTensor = CreateRandomUnitary(Sym, ChargeDimensionsBottom, ChargeDimensionsTop, Legs)
            
            NumberOfLegsBottom = Legs(1);
            NumberOfLegsTop = Legs(2);
            
            OutTensor = SymTensor.CreateTensor(Sym,1,ChargeDimensionsBottom, ChargeDimensionsTop, NumberOfLegsBottom, NumberOfLegsTop);
            
        end
        
        function OutTensor = CreateOnes(Sym, ChargeDimensionsBottom, ChargeDimensionsTop, Legs)
            
            NumberOfLegsBottom = Legs(1);
            NumberOfLegsTop = Legs(2);
            
            OutTensor = SymTensor.CreateTensor(Sym,3,ChargeDimensionsBottom, ChargeDimensionsTop, NumberOfLegsBottom, NumberOfLegsTop);
            
        end
        
        
        function OutTensor = CreateIdentityFrom(InTensor, Legs)
            
            %if ~InTensor.SymHandle.IsSelfDual
                InverseIrrep = InTensor.SymHandle.getInverseIrrep;
                if isequal(InTensor.ChargeSide, [-ones([1,Legs(1)]),ones([1,Legs(2)])])
                    FlagCorrectSides = true;
                elseif isequal(InTensor.ChargeSide, [ones([1,Legs(1)]),-ones([1,Legs(2)])])
                    FlagCorrectSides = false;
                    InTensor.ChargeSide = -InTensor.ChargeSide;
                    InTensor.ChargeDirections = -InTensor.ChargeDirections;
                    
                    %for kk = 1:max(InTensor.Structure(:))
                    %    if 
                    %    InTensor.ChargeLabelsInternal = InverseIrrep(InTensor.ChargeLabelsInternal);
                    %    InTensor.ChargeDirections(self.Structure>0) = -InTensor.ChargeDirections(self.Structure>0);
                    %end
                else
                    error('Error: The Tensor SVD doesn''t work unless all the legs of a side are grouped up at any one time')
                end
            %else
                
            %end
            
            [Matrix,Data] = SymTensor.SymTen2Mat(InTensor, Legs);
            
            for nn = 1:numel(Matrix)
                Matrix{nn} = eye(size(Matrix{nn}));
            end
            
            OutTensor = SymTensor.SymMat2Ten(Matrix,Data);
            
            if ~FlagCorrectSides
                OutTensor.ChargeSide = -OutTensor.ChargeSide;
                OutTensor.ChargeDirections = -OutTensor.ChargeDirections;
            end
            
        end
        
        function OutTensor = CreateRawTensors(TensorEntries, Structure, ChargeDirections, ChargeLabelsExternal, ChargeLabelsInternal,...
                MultiplicitiesInternal, SymHandle,ChargeSide,ChargeSideInt, CorrectDimensions,Braiding,BraidingDirection,KeepZeros)
            
            MinNonZero = 10^-14;
            if nargin<11
                Braiding = [];
            end
            
            if nargin<12
                BraidingDirection = [];
            end
            
            if nargin<10
                CorrectDimensions = false;
                %this means we have already accounted for converting from a
                %matrix into a series of parts
            end
            
            if nargin<13
                KeepZeros = false;
            end
            
            ChargeLegDimensions = zeros([size(SymHandle.getDim,1),size(ChargeLabelsExternal,2)]);
            for kk = 1:size(ChargeLabelsExternal,2)
                [ListExternalCharges,FirstUniqueExternalCharge,AllExternalCharge] = unique(ChargeLabelsExternal(:,kk)');
                for ll = 1:length(ListExternalCharges)
                    ChargeLegDimensions(ListExternalCharges(ll),kk) = size(TensorEntries{FirstUniqueExternalCharge(ll)},kk);
                end
            end
            
            %HERA: need to add checks in here
            
            
            
            %here we correct dimensions as required.
            if CorrectDimensions && SymHandle.IsNonAbelian
                
                SideStructure = zeros(size(Structure));
                [~,IndexExt] = sort(-Structure(Structure<0),'ascend');
                [~,IndexInt] = sort(Structure(Structure>0),'ascend');
                
                SideStructure(Structure<0) = ChargeSide(IndexExt);
                SideStructure(Structure>0) = -ChargeSideInt(IndexInt);
                
                PowerTerms = zeros(size(Structure));
                OutSameSide = SideStructure(1,:)==SideStructure(2,:);
                PowerTerms = PowerTerms + 1*repmat(OutSameSide,[2,1]);
                PowerTerms(:,2:end) = PowerTerms(:,2:end) + (~OutSameSide([1,1],2:end)).*(1-2*([ChargeSideInt;ChargeSideInt] ~= SideStructure(:,2:end) ));
                
                PowerTermsInt = PowerTerms(Structure>0); PowerTermsInt = PowerTermsInt(:)';
                PowerTermsInt = PowerTermsInt(IndexInt) + (-1).^OutSameSide(:,2:end);
                
                PowerTermsExt = PowerTerms(Structure<0); PowerTermsExt = PowerTermsExt(:)';
                PowerTermsExt = PowerTermsExt(IndexExt);
                
                Dim = SymHandle.Dim;
                DimCorrections = prod(reshape(Dim([ChargeLabelsExternal,ChargeLabelsInternal],2),size([ChargeLabelsExternal,ChargeLabelsInternal]))...
                    .^(repmat([PowerTermsExt,PowerTermsInt],[size(ChargeLabelsExternal,1),1])/4),2);
                
                for kk = 1:numel(TensorEntries)
                    TensorEntries{kk} = TensorEntries{kk}*DimCorrections(kk);
                end
            end
            
            if ~KeepZeros
                for kk = numel(TensorEntries):-1:1
                    if max(abs(TensorEntries{kk}(:)))<MinNonZero
                        TensorEntries(kk) = [];
                        MultiplicitiesInternal(kk,:) = [];
                        ChargeLabelsInternal(kk,:) = [];
                        ChargeLabelsExternal(kk,:) = [];
                    end
                end
            end
            
            OutTensor = SymTensor(TensorEntries, Structure, ChargeDirections, MultiplicitiesInternal, ChargeLabelsInternal,...
                ChargeLabelsExternal, Braiding,BraidingDirection,[],ChargeLegDimensions,SymHandle,ChargeSide,ChargeSideInt);
        end
        
        function OutTensor = CreateRawTensorFrom(TensorEntries, InTensor, CorrectDimensions,KeepZeros)
            Structure = InTensor.Structure;
            ChargeDirections = InTensor.ChargeDirections;
            MultiplicitiesInternal = InTensor.MultiplicitiesInternal;
            ChargeLabelsInternal = InTensor.ChargeLabelsInternal;
            ChargeLabelsExternal = InTensor.ChargeLabelsExternal;
            SymHandle = InTensor.SymHandle;
            ChargeSide = InTensor.ChargeSide;
            ChargeSideInt = InTensor.ChargeSideInt;
            Braidings = InTensor.Braidings;
            BraidingDirections = InTensor.BraidingDirections;
            clear InTensor;
            
            if nargin<3||isempty(CorrectDimensions)
                CorrectDimensions = false;
                %this means we have already accounted for converting from a
                %matrix into a series of parts
            end
            if nargin<4 || isempty(KeepZeros)
                KeepZeros = false;
                %this means we have already accounted for converting from a
                %matrix into a series of parts
            end
            
            OutTensor = SymTensor.CreateRawTensors(TensorEntries(:), Structure, ChargeDirections, ChargeLabelsExternal, ChargeLabelsInternal,...
                MultiplicitiesInternal, SymHandle, ChargeSide, ChargeSideInt, CorrectDimensions,Braidings,BraidingDirections,KeepZeros);
        end
        
        function OutTensor = CreateDelta(Sym, ChargeDimensions, Legs)
            if numel(Legs)<2; Legs = [Legs,0]; end;
            
            Location = find(ChargeDimensions~=0, 1,'first');
            Numbers = 1:numel(ChargeDimensions);
            OutTensor = SymTensor.CreateOnes(Sym, ChargeDimensions.*(Numbers==Location), ChargeDimensions.*(Numbers==Location), Legs);
            Temp = zeros(size(OutTensor.TensorEntries{1}));
            mm = sum((repmat(1:size(Temp,1), [numel(size(Temp)),1])-1).*repmat(size(Temp,1).^((0:(numel(size(Temp))-1))'), [1,size(Temp,1)]),1)+1;
            Temp(mm) = ones([1,size(Temp,1)]);
            
            for kk = 1:numel(OutTensor.TensorEntries)
                OutTensor.TensorEntries{kk} = Temp;
            end
            Location = Location+1;
            
            while Location<=numel(ChargeDimensions) 
                
                if ChargeDimension(Location) == 0
                    Location = Location+1;
                    continue;
                end
                Location = Location+1;
                
                TempTensor = SymTensor.CreateOnes(Sym, ChargeDimensions.*(Numbers==Location), ChargeDimensions.*(Numbers==Location), Legs);
                Temp = zeros(size(TempTensor.TensorEntries{1}));
                mm = (repmat(1:size(Temp,1), [1,size(Temp,1)])-1).*repmat(size(Temp,1).^(0:(size(Temp,1)-1)), [1,size(Temp,1)])+1;
                Temp(mm) = ones([1,size(Temp,1)]);
                 
                for kk = 1:numel(TempTensor.TensorEntries)
                    TempTensor.TensorEntries{kk} = Temp;
                end
                
                OutTensor = OutTensor+TempTensor;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Load Ratio
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %this signals the updating with charges
        
        function [Matrix, Data] = SymTen2Mat(Tensor, Legs, Shape)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first step is to work out how we want to reshape it.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ReshapeData = 0;
            
            if nargin<2 | isempty(Legs)
                NumberLegsBottom = sum(Tensor.ChargeSide ==-1);
                NumberLegsTop = sum(Tensor.ChargeSide ==+1);
            else
                NumberLegsBottom = Legs(1);
                NumberLegsTop = Legs(2);
            end
            
            
            if nargin==3
                
                if SymCanInsert(Tensor,Shape)
                    ChargeLabelsA = [Tensor.ChargeLabelsExternal,Tensor.ChargeLabelsInternal,Tensor.MultiplicitiesInternal];
                    ChargeLabelsB = [Shape.ChargeLabelsExternal,Shape.ChargeLabelsInternal,Shape.MultiplicitiesInternal];
                    
                    [~,Index,FullIndex] = unique([ChargeLabelsB;ChargeLabelsA],'rows','first');
                    
                    if any(Index>size(ChargeLabelsB,1))
                        error('Affirmation Error: Something is wrong with Shape and working out unique charges')
                    end
                    
                    Locations = Index(FullIndex((size(ChargeLabelsB,1)+1):end));
                    
                    ReplaceTensorEntries = Shape.TensorEntries;
                    if Shape.FlagIsAShape
                        for kk = 1:numel(ReplaceTensorEntries)
                           ReplaceTensorEntries{kk} = zeros(Shape.TensorEntriesSizes(kk,:));
                        end
                    else
                        for kk = 1:numel(ReplaceTensorEntries)
                           ReplaceTensorEntries{kk} = ReplaceTensorEntries{kk}*0;
                        end
                    end
                    
                    ReplaceTensorEntries(Locations) = Tensor.TensorEntries;
                    
                    Tensor = SymTensor.CreateRawTensorFrom(ReplaceTensorEntries, Shape,[],true);
                    
                else
                    error('Error: The input Shape is either not larger then the tensor we wish to find the matrix of, or of a qualitatively different layout, therefore we cannot insert Tensor Entries into Shape');
                end
            end
            
            %check that this agrees with my code:
            
            %check we haven't done something stupid, only way we could is
            %by putting the wrong number of legs in:
            if (NumberLegsBottom + NumberLegsTop) ~= size(Tensor.ChargeLabelsExternal,2)
                error('The number of input legs is wrong for this tensor, we should have a different number')
            end
            
            if sum(Tensor.ChargeSide ==+1) ~= NumberLegsTop;
                error('Error: Wrong number of Top Legs')
            end
            
            if sum(Tensor.ChargeSide ==-1) ~= NumberLegsBottom;
                error('Error: Wrong number of Bottom Legs')
            end
            
            if ~isequal(Tensor.ChargeSide, [ones([1,NumberLegsBottom]),-ones([1,NumberLegsTop])])
                [~,PermuteOrder] = sort(Tensor.ChargeSide);
                Tensor = Tensor.Permute(PermuteOrder);
            end
            
            
%            Tensor = Tensor.SymReshape(ReshapeData);
            
            Tensor = Tensor.SortLabels;
            
            
            %check that the shape is correct
            NumbersLegs = (1:(NumberLegsTop+NumberLegsBottom));
            
            if NumberLegsBottom~=1
                if ~isequal(sort(Tensor.StoredLocations{Tensor.Structure(1,1)},'ascend'), NumbersLegs(Tensor.ChargeSide == -1))
                    error('Error: The organisation sets the Bottom to be wrong')
                end
            else
                if Tensor.Structure(1,1) ~= -1;
                    error('Error: The organisation sets the Bottom to be wrong')
                end
            end
            
            if NumberLegsTop~=1
                if ~isequal(sort(Tensor.StoredLocations{Tensor.Structure(2,1)},'ascend'), NumbersLegs(Tensor.ChargeSide == +1))
                    error('Error: The organisation sets the Top to be wrong')
                end
            else
                if Tensor.Structure(2,1) ~= -(1+NumberLegsBottom);
                    error('Error: The organisation sets the Top to be wrong')
                end
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %second step is to sort into matrix blocks
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                Dim = Tensor.SymHandle.Dim;
            %{
            if Tensor.SymHandle.IsNonAbelian
                
                SideStructure = zeros(size(Tensor.Structure));
                [~,IndexExt] = sort(-Tensor.Structure(Tensor.Structure<0),'ascend');
                [~,IndexInt] = sort(Tensor.Structure(Tensor.Structure>0),'ascend');
                
                SideStructure(Tensor.Structure<0) = Tensor.ChargeSide(IndexExt);
                SideStructure(Tensor.Structure>0) = -Tensor.ChargeSideInt(IndexInt);
                
                PowerTerms = zeros(size(Tensor.Structure));
                OutSameSide = SideStructure(1,:)==SideStructure(2,:);
                PowerTerms = PowerTerms + 1*repmat(OutSameSide,[2,1]);
                PowerTerms(:,2:end) = PowerTerms(:,2:end) + (~OutSameSide([1,1],2:end)).*(1-2*([Tensor.ChargeSideInt;Tensor.ChargeSideInt] ~= SideStructure(:,2:end) ));
                
                PowerTermsInt = PowerTerms(Tensor.Structure>0); PowerTermsInt = PowerTermsInt(:)';
                PowerTermsInt = PowerTermsInt(IndexInt) + (-1).^OutSameSide(:,2:end);
                
                PowerTermsExt = PowerTerms(Tensor.Structure<0); PowerTermsExt = PowerTermsExt(:)';
                PowerTermsExt = PowerTermsExt(IndexExt);
                
                Dim = Tensor.SymHandle.Dim;
                DimCorrections = prod(reshape(Dim([Tensor.ChargeLabelsExternal,Tensor.ChargeLabelsInternal],2),size([Tensor.ChargeLabelsExternal,Tensor.ChargeLabelsInternal]))...
                    .^(-repmat([PowerTermsExt,PowerTermsInt],[size(Tensor.ChargeLabelsExternal,1),1])/4),2);
                
                for kk = 1:numel(Tensor.TensorEntries)
                    Tensor.TensorEntries{kk} = Tensor.TensorEntries{kk}*DimCorrections(kk);
                end
            end
            %}
            
            if NumberLegsBottom*NumberLegsTop ~= 0
                
                if NumberLegsBottom > 1
                    [UniqueEntries,~,BlockLabelBottom] = unique(Tensor.ChargeLabelsInternal(:,NumberLegsBottom-1));
                else %NumberLegsBottom == 1; as it can't be zero
                    [UniqueEntries,~,BlockLabelBottom] = unique(Tensor.ChargeLabelsExternal(:,1));
                end
                
                
                
                %UniqueEntries is the list of unique values that the matrix
                %can be valued for (the rep when in matrix format), we know
                %that we have all of them because we have to have the
                %structure consistant so if they exist on the bottom they
                %must also exist on the top (B<T), if they exist on the top
                %then they must have a rep on the bottom and so we still
                %pick it up (T<B) So we get all of the reps we need to care
                %about in this
                
                %don't forget we are in the form of a tree from the centre
                %with the first Legs(1) are the bottom, the next Legs(2)
                %are the top legs
                
                
                if NumberLegsTop > 1
                    [UniqueEntries2,~,BlockLabelTop] = unique(Tensor.ChargeLabelsInternal(:,end));
                else %NumberLegsTop == 1; as it can't be zero
                    [UniqueEntries2,~,BlockLabelTop] = unique(Tensor.ChargeLabelsExternal(:,end));
                end
                
                ChargeLegDimTemp = Tensor.ChargeLegDimensions;
                
                
                
                
                
                for jj = 1:size(UniqueEntries,1)
                    
                    KeptLabelsBottom = BlockLabelBottom == jj;
                    
                    
                    if NumberLegsBottom > 1
                        
                        BottomLabels{jj} = [Tensor.ChargeLabelsExternal(KeptLabelsBottom,1:(NumberLegsBottom)),...
                                    Tensor.ChargeLabelsInternal(KeptLabelsBottom,1:(NumberLegsBottom-1)),...
                                    Tensor.MultiplicitiesInternal(KeptLabelsBottom,1:(NumberLegsBottom-1))];
                        
                        
                        [UniqueEntriesBottom{jj},~,LabelsBottom{jj}] = unique(BottomLabels{jj},'rows');
                        TT = Dim(:,2);
                        FactorBottom{jj} = (prod(TT(UniqueEntriesBottom{jj}(:,1:NumberLegsBottom)),2)).^(1/4);
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Load Ratio
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                    else %NumberLegsBottom == 1; as it can't be zero
                        
                        %now we know that it has to be the Label,
                        
                        UniqueEntriesBottom{jj} = UniqueEntries(jj);
                        
                        LabelsBottom{jj} = ones([sum(KeptLabelsBottom),1]);
                        
                        FactorBottom{jj} = (Dim(UniqueEntries(jj),2)).^(1/4);
                    end
                    
                    
                    DimensionsBottom{jj} = zeros([size(UniqueEntriesBottom{jj},1), NumberLegsBottom]);
                    
                    for kk = 1:size(UniqueEntriesBottom{jj},1)
                        for ll = 1:NumberLegsBottom
                            
                            DimensionsBottom{jj}(kk,ll) = ChargeLegDimTemp(UniqueEntriesBottom{jj}(kk,ll),ll);
                            
                        end
                    end
                    
                    DimensionsBottomSingle{jj} = prod(DimensionsBottom{jj}, 2);
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    KeptLabelsTop = BlockLabelTop == jj;
                    
                    
                    if NumberLegsTop > 1
                        
                        TopLabels{jj} = [Tensor.ChargeLabelsExternal(KeptLabelsTop,NumberLegsBottom+(1:NumberLegsTop)),...
                                    Tensor.ChargeLabelsInternal(KeptLabelsTop,NumberLegsBottom-1+(1:(NumberLegsTop-1))),...
                                    Tensor.MultiplicitiesInternal(KeptLabelsTop,NumberLegsBottom-1+(1:(NumberLegsTop-1)))]; %extra plus one in the middle to account for the useless Multiplicitiy for the zero source
                        
                        
                        [UniqueEntriesTop{jj},~,LabelsTop{jj}] = unique(TopLabels{jj},'rows');
                        TT = Dim(:,2);
                        FactorTop{jj} = (prod(TT(UniqueEntriesTop{jj}(:,1:NumberLegsTop)),2)).^(1/4);
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Load Ratio
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                    else %NumberLegsTop == 1; as it can't be zero
                        
                        UniqueEntriesTop{jj} = UniqueEntries(jj);
                        LabelsTop{jj} = ones([sum(KeptLabelsTop),1]);
                        
                        FactorTop{jj} = (Dim(UniqueEntries(jj),2)).^(1/4);
                    end
                    
                    DimensionsTop{jj} = zeros([size(UniqueEntriesTop{jj},1), NumberLegsTop]);
                    
                    for kk = 1:size(UniqueEntriesTop{jj},1)
                        for ll = 1:NumberLegsTop
                            
                            DimensionsTop{jj}(kk,ll) = ChargeLegDimTemp(UniqueEntriesTop{jj}(kk,ll),ll+NumberLegsBottom);
                            
                        end
                    end
                    
                    DimensionsTopSingle{jj} = prod(DimensionsTop{jj}, 2);
                    
                end
                
            else
                %then one of them is zero which means that we only have one
                %block of trivial charges
                
                BlockLabelBottom = ones([size(Tensor.ChargeLabelsInternal,1),1]);
                BlockLabelTop = ones([size(Tensor.ChargeLabelsInternal,1),1]);
                
                UniqueEntries = Tensor.SymHandle.GetTrivialIrrep;
                %if either of them are zero then the only possible charges
                %are the trivial irreps.
                
                
                if NumberLegsTop == 0
                    DimensionsTop = 1;
                    DimensionsTopSingle = 1;
                    FlagOnlyLoop = ones([size(Tensor.ChargeLabelsExternal,1),1]);
                    
                    BottomLabels = [Tensor.ChargeLabelsExternal, Tensor.ChargeLabelsInternal, Tensor.MultiplicitiesInternal];
                    TopLabels = zeros([1,0]);
                    
                    [UniqueEntriesBottom,~,LabelsBottom] = unique(BottomLabels,'rows');
                    
                    DimensionsBottom = zeros([size(UniqueEntriesBottom,1), NumberLegsBottom]);
                    
                    ChargeLegDimTemp = Tensor.ChargeLegDimension;
                    
                    for kk = 1:UniqueEntriesBottom
                        for ll = 1:NumberLegsBottom
                            
                            DimensionsBottom(kk,ll) = ChargeLegDimTemp(UniqueEntriesBottom(kk,ll),ll);
                            
                        end
                    end
                    
                    DimensionsBottomSingle = prod(DimensionBottom, 2);
                    
                    
                    
                else %NumberLegsBottom == 0
                    DimensionsBottom = 1;
                    DimensionsBottomSingle = 1;
                    LabelsBottom = ones([size(Tensor.ChargeLabelsExternal,1),1]);
                    
                    TopLabels = [Tensor.ChargeLabelsExternal, Tensor.ChargeLabelsInternal, Tensor.MultiplicitiesInternal];
                    BottomLabels = zeros([1,0]);
                    [UniqueEntriesTop,~,LabelsTop] = unique(TopLabels,'rows');
                    
                    DimensionsTop = zeros([size(UniqueEntriesTop,1), NumberLegsTop]);
                    
                    ChargeLegDimTemp = Tensor.ChargeLegDimension;
                    
                    for kk = 1:UniqueEntriesTop
                        for ll = 1:NumberLegsTop
                            
                            DimensionsTop(kk,ll) = ChargeLegDimTemp(UniqueEntriesTop(kk,ll),ll+NumberLegsBottom);
                            
                        end
                    end
                    
                    DimensionsTopSingle = prod(DimensionTop, 2);
                    
                end
                
                %now we need to make this one copy
                DimensionsTop = {DimensionsTop};
                DimensionsTopSingle = {DimensionsTopSingle};
                UniqueEntriesTop = {UniqueEntriesTop};
                DimensionsBottom = {DimensionsBottom};
                DimensionsBottomSingle = {DimensionsBottomSingle};
                UniqueEntriesBottom = {UniqueEntriesBottom};
                
                
            end
            
            
            

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Third step is to sort blocks and then put together
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            
            
            %make matricies, note first we need to fill if empty
            
            MatrixParts = cell([size(UniqueEntries,1),1]);
                
            for jj = 1:size(UniqueEntries,1)
                
                MatrixParts{jj} = cell([max(LabelsBottom{jj},[],1),max(LabelsTop{jj},[],1)]);
                
                TensorParts = Tensor.TensorEntries(BlockLabelBottom==jj);
                
                Counter = 1;
                
                for kk = [LabelsBottom{jj},LabelsTop{jj}]'%[UniqueEntriesBottom{jj},UniqueEntriesTop{jj}]'
                    MatrixParts{jj}{kk(1),kk(2)} = reshape(TensorParts{Counter}, [DimensionsBottomSingle{jj}(kk(1)),DimensionsTopSingle{jj}(kk(2))])...
                        /sqrt(Dim(UniqueEntries(jj),2))*...
                        ((prod([Dim(UniqueEntriesBottom{jj}(kk(1),1:NumberLegsBottom),2);Dim(UniqueEntriesTop{jj}(kk(2),1:NumberLegsTop),2)]))^(1/4));
                        
                %/sqrt(Dim(UniqueEntries(jj),2))*(FactorBottom{jj}(kk(1))*FactorTop{jj}(kk(2)));
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Load Ratio Ten2Mat
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    Counter = Counter+1;
                end
                
                
                
                EmptyEntries = find(  (cellfun('isempty', MatrixParts{jj}(:))));
                
                kk1 = mod(EmptyEntries-1,size(MatrixParts{jj},1))+1;
                kk2 = round((EmptyEntries-kk1)/size(MatrixParts{jj},1));
                kk1 = kk1(:); kk2 = kk2(:)+1;
                
                for kk = [kk1,kk2]'
                    MatrixParts{jj}{kk(1),kk(2)} = zeros([DimensionsBottomSingle{jj}(kk(1)),DimensionsTopSingle{jj}(kk(2))]);
                end
                
                Matrix{jj} = cell2mat(MatrixParts{jj});
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %forth step is to work out what data we need
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Data = {UniqueEntries, DimensionsBottom, DimensionsTop, Legs, UniqueEntriesBottom,...
                UniqueEntriesTop, Tensor.ChargeDirections, Tensor.Braidings, Tensor.BraidingDirections,...
                Tensor.SymHandle,Tensor.Structure, ReshapeData};
            
            
            
        end
        
        function Tensor = SymMat2Ten(Matrix,Data,KeepZeros)
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %zero step is to check that nothing stupid is being done
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if nargin<3
                KeepZeros = false;
            end
            
            AllowedError = 10^-15;
            if KeepZeros
                AllowedError = -1;
            end
            
            if ~isa(Data,'cell')
                error('Error: Data should be a cell of 12 entries')
            end
            
            if numel(Data)<12
                error('Error: Data should be a cell of at least 12 entries')
            end
            
            
            UniqueEntries = Data{1};
            
            DimensionsBottom = Data{2};
            DimensionsTop = Data{3};
            DimensionsBottomSingle = cell(size(DimensionsBottom));
            DimensionsTopSingle = cell(size(DimensionsTop));
            for jj = 1:numel(DimensionsBottom)
                DimensionsBottomSingle{jj} = prod(DimensionsBottom{jj},2);
                DimensionsTopSingle{jj} = prod(DimensionsTop{jj},2);
                
                if sum(DimensionsBottomSingle{jj},1) ~= size(Matrix{jj},1)
                    error('Error: Dimensions Bottom (entry 2) doesn''t agree with the size of the matrix input')
                end
                if sum(DimensionsTopSingle{jj},1) ~= size(Matrix{jj},2)
                    error('Error: Dimensions Top (entry 3) doesn''t agree with the size of the matrix input')
                end
                
            end
            Legs = Data{4};
            
            UniqueEntriesBottom = Data{5};
            UniqueEntriesTop = Data{6};
            
            ChargeDirections = Data{7};
            Braidings = Data{8};
            BraidingDirections = Data{9};
            SymHandle = Data{10};
            Structure = Data{11};
            

            
            ReshapeData = Data{12};
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Zeroth step is to introduce checks 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %check the symmetry first
            if ~isa(SymHandle, 'Symmetry')
                error('Error: Data is wrong. This is not a Symmetry (entry 10)')
            end
            
            TempDim = SymHandle.getDimExternal;
            Dim = SymHandle.getDim;
            
            
            %now check Legs
            if ~IsInteger(Legs)
                error('Error: Data is wrong. The Legs (entry 4) should be integers.')
            end
            
            if numel(DimensionsBottom)~=numel(UniqueEntries)
                error('Error: Data is wrong. The Dimensions Bottom (entry 2) number of entries does not agree with the number of entries in UniqueEntries (entry 1)')
            end
            for jj = 1:numel(DimensionsBottom)
                if ~IsInteger(DimensionsBottom{jj})
                    error('Error: Data is wrong. The Dimensions Bottom (entry 2) should be integers.')
                end
                if Legs(1)~=size(DimensionsBottom{jj},2)
                    error('Error: Data is wrong. The Legs (entry 4) don''t agree with the Dimensions Bottom (entry 2)')
                end
                if length(size(DimensionsBottom{jj})) > 2
                    error('Error: Data is wrong. The Dimensions Bottom(entry 2) has too many dimensions')
                end
            end
            
            if numel(DimensionsTop)~=numel(UniqueEntries)
                error('Error: Data is wrong. The Dimensions Top (entry 3) number of entries does not agree with the number of entries in UniqueEntries (entry 1)')
            end
            for jj = 1:numel(DimensionsTop)
                if ~IsInteger(DimensionsTop{jj})
                    error('Error: Data is wrong. The Dimensions Top (entry 3) should be integers.')
                end
                if Legs(2)~=size(DimensionsTop{jj},2)
                    error('Error: Data is wrong. The Legs (entry 4) don''t agree with the Dimensions Top (entry 3)')
                end
                if length(size(DimensionsTop{jj})) > 2
                    error('Error: Data is wrong. The Dimensions Top(entry 3) has too many dimensions')
                end
            end
            
            LegsInternal = [0,0];
            
            if Legs(1) < 2 %so 0 or 1, when there is no internal
                AllNumbersBottom = Legs(1);
                LegsInternal(1) = 0;
            else
                AllNumbersBottom = 3*Legs(1)-2;
                LegsInternal(1) = Legs(1)-1;
            end
            
            
            if numel(UniqueEntriesBottom)~=numel(UniqueEntries)
                error('Error: Data is wrong. The Unique Entries Bottom (entry 5) number of entries does not agree with the number of entries in UniqueEntries (entry 1)')
            end
            for jj = 1:numel(UniqueEntriesBottom)
                if ~IsInteger(UniqueEntriesBottom{jj})
                    error('Error: Data is wrong. The Unique Entries Bottom (entry 5) should be integers.')
                end
                if AllNumbersBottom~=size(UniqueEntriesBottom{jj},2)
                    error('Error: Data is wrong. The Legs (entry 4) don''t agree with the Unique Entries Bottom (entry 5)')
                end
                if length(size(UniqueEntriesBottom{jj})) > 2
                    error('Error: Data is wrong. The Unique Entries Bottom (entry 5) has too many dimensions')
                end
            end
            
            
            
            if Legs(2) < 2 %so 0 or 1, when there is no internal
                AllNumbersTop = Legs(2);
                LegsInternal(2) = 0;
            else
                AllNumbersTop = 3*Legs(2)-2;
                LegsInternal(2) = Legs(2)-1;
            end
            
            if numel(UniqueEntriesTop)~=numel(UniqueEntries)
                error('Error: Data is wrong. The Unique Entries Top (entry 6) number of entries does not agree with the number of entries in UniqueEntries (entry 1)')
            end
            for jj = 1:numel(UniqueEntriesTop)
                if ~IsInteger(UniqueEntriesTop{jj})
                    error('Error: Data is wrong. The Unique Entries Top (entry 6) should be integers.')
                end
                if AllNumbersTop~=size(UniqueEntriesTop{jj},2)
                    error('Error: Data is wrong. The Legs (entry 4) don''t agree with the UniqueEntries Top (entry 6)')
                end
                if length(size(UniqueEntriesBottom{jj})) > 2
                    error('Error: Data is wrong. The Unique Entries Bottom (entry 6) has too many dimensions')
                end
            end
            
            if Legs(1)+Legs(2)==0
                if ~isempty(Structure)
                    error('Error: Structure (entry 11) is wrong, it should be empty for a zero entry number')
                end
                if ~isempty(ChargeDirections)
                    error('Error: Charge Directions (entry 7) is wrong, it should be empty for a zero entry number')
                end
            else
                if ~IsInteger(Structure)
                    error('Error: Data is wrong. The Structure (entry 11) should be integers.')
                end
                if ~IsInteger(ChargeDirections)
                    error('Error: Data is wrong. The Charge Directions (entry 7) should be integers.')
                end
                
                if Legs(1)+Legs(2) == 2
                    StructureSorted = [-2;-1];
                else
                    StructureSorted = [-(Legs(1)+Legs(2)):-1, 1:(Legs(1)+Legs(2)-2)]';
                end
                
                if ~isequal(sort(Structure(:),'ascend'), StructureSorted) || ~isequal(size(Structure), [2,Legs(1)+Legs(2)-1])
                    error('Error: Data is wrong. The structure (entry 11) is wrong')
                end
                
                if any(ChargeDirections(:)~=1&ChargeDirections(:)~=-1) || ~isequal(size(ChargeDirections), [2,Legs(1)+Legs(2)-1])
                    error('Error: Data is wrong. The Charge Directions (entry 7) is wrong')
                end
                
            end
            
            
            %now check UniqueEntriesBottom and UniqueEntriesTop to make
            %sure they agree with the Symmetry
            
            %WILL IMPLIMENT LATER
            
            
            %now check UniqueEntries to make sure they agree with the
            %Symmetry and the UniqueEntriesBottom and UniqueEntriesTop
            
            %WILL IMPLIMENT LATER
            
            
            %check that the braiding details are sensible
            
            %WILL IMPLIMENT LATER
            
            
            %check that the reshape data makes sense
            
            %WILL IMPLIMENT LATER
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %first step is to work out how we want to reshape it.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ReshapeData = 0;
            InverseOfReshapeData = 0;
            
            NumberLegsBottom = Legs(1);
            NumberLegsTop = Legs(2);

            %DimensionsBottom{jj} is from UniqueEntriesBottom
            TempDim = SymHandle.getDimExternal;
            ChargeLegDimensions = zeros([size(TempDim,1), NumberLegsBottom+NumberLegsTop]);
            
            for jj = 1:numel(UniqueEntriesTop)
                %go through each block of the matrix
                for kk1 = (TempDim(:,1))'
                    for kk2 = 1:NumberLegsBottom
                        
                        Temp = find(UniqueEntriesBottom{jj}(:,kk2) == kk1);
                        if ~isempty(Temp)
                            Dimension = DimensionsBottom{jj}(Temp,kk2);
                            
                            if all(Dimension(:)==repmat(Dimension(1),[numel(Dimension),1]))
                                %then we are safe
                                if ChargeLegDimensions(kk1,kk2) == 0
                                    ChargeLegDimensions(kk1,kk2) = Dimension(1);
                                else
                                    if ChargeLegDimensions(kk1,kk2) ~= Dimension(1)
                                        error('Error: We should have all the legs agreeing to what the Dimension should be in DimensionBottom (entry 2)')
                                    end
                                end
                            else
                                error('Error: We should have all the legs agreeing to what the Dimension should be in DimensionBottom (entry 2)')
                            end
                        end
                        
                    end
                    
                    for kk2 = 1:NumberLegsTop
                        Temp = find(UniqueEntriesTop{jj}(:,kk2) == kk1);
                        
                        if ~isempty(Temp)
                            Dimension = DimensionsTop{jj}(Temp,kk2);
                            
                            if all(Dimension(:)==repmat(Dimension(1),[numel(Dimension),1]))
                                %then we are safe
                                if ChargeLegDimensions(kk1,kk2+NumberLegsBottom) == 0
                                    ChargeLegDimensions(kk1,kk2+NumberLegsBottom) = Dimension(1);
                                else
                                    if ChargeLegDimensions(kk1,kk2+NumberLegsBottom) ~= Dimension(1)
                                        error('Error: We should have all the legs agreeing to what the Dimension should be in DimensionTop (entry 3)')
                                    end
                                end
                            else
                                error('Error: We should have all the legs agreeing to what the Dimension should be in DimensionTop (entry 3)')
                            end
                        end
                    end
                end
            end
            
            
            %Now generate the Unique components of external and internal
            %charges and multiplicities
            
            Fusion2 = SymHandle.getFusion2;
            TrivialIrrep = SymHandle.getTrivialIrrep;
            InverseIrrep = SymHandle.getInverseIrrep;
            
            
            
            
            
            
            %here we work out exactly what we should be expecting to see
            %for the entry in UniqueEntriesBottom(Top) and UniqueEntries
            %i.e. if the charge is going in the opposite direction to the
            %expected one then we have to inverse it so that we can check.
            if NumberLegsBottom*NumberLegsTop ~=0
                
                if ChargeDirections(1,1) == -1
                    UniqueEntriesVeryBottom = UniqueEntries;
                else
                    UniqueEntriesVeryBottom = InverseIrrep(UniqueEntries);
                end
                
                if ChargeDirections(2,1) == 1
                    UniqueEntriesVeryTop = UniqueEntries;
                else
                    UniqueEntriesVeryTop = InverseIrrep(UniqueEntries);
                end
                
            else
                UniqueEntriesVeryBottom = UniqueEntries;
                UniqueEntriesVeryTop = UniqueEntries;
            end
            
            
            %Now break down the input entries into the
            %ChargeLabelsExternal, ChargeLabelsInternal and
            %MultiplicitiesInternal
            for jj = 1:numel(UniqueEntries);
                
                if NumberLegsBottom>1
                    ChargeLabelsExternalUniqueBottom{jj} = UniqueEntriesBottom{jj}(:,1:NumberLegsBottom);
                    ChargeLabelsInternalUniqueBottom{jj} = UniqueEntriesBottom{jj}(:,(1:(NumberLegsBottom-1))+NumberLegsBottom);
                    MultiplicitiesInternalUniqueBottom{jj} = UniqueEntriesBottom{jj}(:,(1:(NumberLegsBottom-1))+2*NumberLegsBottom -1);
                    
                    if any(ChargeLabelsInternalUniqueBottom{jj}(:,end)~=UniqueEntriesVeryBottom(jj))
                        error('The charges at the bottom (entry 5) should agree with the UniqueEntries (entry 1)')
                    end
                elseif NumberLegsBottom == 1
                    ChargeLabelsExternalUniqueBottom{jj} = UniqueEntriesBottom{jj};
                    ChargeLabelsInternalUniqueBottom{jj} = zeros([1,0]);
                    MultiplicitiesInternalUniqueBottom{jj} = zeros([1,0]);
                    
                    if any(ChargeLabelsExternalUniqueBottom{jj}(:,end)~=UniqueEntriesVeryBottom(jj))
                        error('The charges at the bottom (entry 5) should agree with the UniqueEntries (entry 1)')
                    end
                else %if NumberLegsBottom == 0
                    ChargeLabelsExternalUniqueBottom{jj} = zeros([1,0]);
                    ChargeLabelsInternalUniqueBottom{jj} = zeros([1,0]);
                    MultiplicitiesInternalUniqueBottom{jj} = zeros([1,0]);
                    if jj>1
                        error('There are no bottom degrees of freedom, this should only give a trival bottom')
                    elseif UniqueEntriesVeryTop ~= TrivialIrrep
                        error('There are no bottom degrees of freedom, this should only give a trival bottom')
                    end
                end
                
                if NumberLegsTop>1
                    ChargeLabelsExternalUniqueTop{jj} = UniqueEntriesTop{jj}(:,1:NumberLegsTop);
                    ChargeLabelsInternalUniqueTop{jj} = UniqueEntriesTop{jj}(:,(1:(NumberLegsTop-1))+NumberLegsTop);
                    MultiplicitiesInternalUniqueTop{jj} = UniqueEntriesTop{jj}(:,(1:(NumberLegsTop-1))+2*NumberLegsTop -1);
                    
                    if any(ChargeLabelsInternalUniqueTop{jj}(:,end)~=UniqueEntriesVeryTop(jj))
                        error('The charges at the top (entry 6) should agree with the UniqueEntries (entry 1)')
                    end
                    
                elseif NumberLegsTop == 1
                    ChargeLabelsExternalUniqueTop{jj} = UniqueEntriesTop{jj};
                    ChargeLabelsInternalUniqueTop{jj} = zeros([1,0]);
                    MultiplicitiesInternalUniqueTop{jj} = zeros([1,0]);
                    
                    if any(ChargeLabelsExternalUniqueTop{jj}(:,end)~=UniqueEntriesVeryTop(jj))
                        error('The charges at the top (entry 6) should agree with the UniqueEntries (entry 1)')
                    end
                else %ifNumberLegsTop == 0
                    ChargeLabelsExternalUniqueTop{jj} = zeros([1,0]);
                    ChargeLabelsInternalUniqueTop{jj} = zeros([1,0]);
                    MultiplicitiesInternalUniqueTop{jj} = zeros([1,0]);
                    
                    if jj>1
                        error('There are no top degrees of freedom, this should only give a trival top')
                    elseif UniqueEntriesVeryTop ~= TrivialIrrep
                        error('There are no top degrees of freedom, this should only give a trival top')
                    end
                                       
                    
                end
                
                
                %verify
                
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %second step is to convert the matrix into matrix parts and
            %drop if zeros
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            MatrixParts = cell([size(Matrix,1),1]);
            
            TempTensorEntries = cell([0,1]);
            TempTensorLabels = zeros([0,2]);
            TempChargeLabelsExternal = [];
            TempChargeLabelsInternal = [];
            TempMultiplicitiesInternal = [];
            
            for jj = 1:length(Matrix)
                
                
                MatrixParts{jj} = mat2cell(Matrix{jj},DimensionsBottomSingle{jj},DimensionsTopSingle{jj});
                
                MatrixPartsNumbers = [repmat(1:size(MatrixParts{jj},1),[1,size(MatrixParts{jj},2)]);...
                    reshape(repmat((1:size(MatrixParts{jj},2)),[size(MatrixParts{jj},1),1]),[1,size(MatrixParts{jj},1)*size(MatrixParts{jj},2)])];
                
               
                for kk = MatrixPartsNumbers
                    if any(any(abs(MatrixParts{jj}{kk(1),kk(2)})>AllowedError,1),2) %if non-zero
                        MatrixParts{jj}{kk(1),kk(2)} = reshape(MatrixParts{jj}{kk(1),kk(2)}, [DimensionsBottomSingle{jj}(kk(1),:),DimensionsTopSingle{jj}(kk(2),:)])...
                            *sqrt(Dim(UniqueEntries(jj),2))/...
                            ((prod([Dim(UniqueEntriesBottom{jj}(kk(1),1:NumberLegsBottom),2);Dim(UniqueEntriesTop{jj}(kk(2),1:NumberLegsTop),2)]))^(1/4));
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Load Ratio Mat2Ten
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    else
                        MatrixParts{jj}{kk(1),kk(2)} = [];
                    end
                end
            
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Third step is to convert the input into corresponding
            %ChargeLabels and multiplicities
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                FullEntries = find(not(cellfun('isempty', MatrixParts{jj}(:))));
                
                CurrentLength = size(TempTensorEntries,1);
                
                if length(FullEntries)>0
                    TempTensorEntries{CurrentLength+length(FullEntries),1} = [];
                    
                    kk1 = mod(FullEntries-1,size(MatrixParts{jj},1))+1;
                    kk2 = round((FullEntries-kk1)/size(MatrixParts{jj},1));
                    kk1 = kk1(:); kk2 = kk2(:)+1;
                    TempTensorEntries(CurrentLength + (1:length(FullEntries)),1) = reshape(MatrixParts{jj}(FullEntries),[length(FullEntries),1]);
                    
                    %fix up external charges
                    if ~isempty(TempChargeLabelsExternal)
                        TempChargeLabelsExternal = [TempChargeLabelsExternal;...
                            [ChargeLabelsExternalUniqueBottom{jj}(kk1,:),ChargeLabelsExternalUniqueTop{jj}(kk2,:)]];
                    else
                        TempChargeLabelsExternal = [ChargeLabelsExternalUniqueBottom{jj}(kk1,:),ChargeLabelsExternalUniqueTop{jj}(kk2,:)];
                    end
                    
                    
                    %fix up internal charges
                    if ~isempty(TempChargeLabelsInternal)
                        TempChargeLabelsInternal = [TempChargeLabelsInternal;...
                        [ChargeLabelsInternalUniqueBottom{jj}(kk1,:),ChargeLabelsInternalUniqueTop{jj}(kk2,:)]];
                    else
                        TempChargeLabelsInternal = [ChargeLabelsInternalUniqueBottom{jj}(kk1,:),ChargeLabelsInternalUniqueTop{jj}(kk2,:)];
                    end
                    
                    %fix up Multiplicity charges
                    if NumberLegsBottom*NumberLegsTop ~= 0 %if we need to worry about the useless middle Multiplicity of ones
                        if ~isempty(TempMultiplicitiesInternal)
                            TempMultiplicitiesInternal = [TempMultiplicitiesInternal;...
                                [MultiplicitiesInternalUniqueBottom{jj}(kk1,:),MultiplicitiesInternalUniqueTop{jj}(kk2,:)]];
                        else
                            TempMultiplicitiesInternal = [MultiplicitiesInternalUniqueBottom{jj}(kk1,:),MultiplicitiesInternalUniqueTop{jj}(kk2,:)];
                        end
                    else %one side has no legs out.
                        if ~isempty(TempMultiplicitiesInternal)
                            TempMultiplicitiesInternal = [TempMultiplicitiesInternal;...
                                [MultiplicitiesInternalUniqueBottom{jj}(kk1,:),MultiplicitiesInternalUniqueTop{jj}(kk2,:)]];
                        else
                            TempMultiplicitiesInternal = [MultiplicitiesInternalUniqueBottom{jj}(kk1,:),MultiplicitiesInternalUniqueTop{jj}(kk2,:)];
                        end
                    end
                end
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            for ll = 1:size(TempTensorEntries,1)
                if ~isempty(TempTensorEntries{ll})
                   Dimensions = ChargeLegDimensions(TempChargeLabelsExternal(ll,:)+((1:size(TempChargeLabelsExternal,2))-1)*size(ChargeLegDimensions,1));
                   TempTensorEntries{ll} = reshape(TempTensorEntries{ll},Dimensions(:)');
                end
            end
            
            
            
            ChargeSide = [-ones([1,Legs(1)]),ones([1,Legs(2)])];
            
            ChargeSideInt = [-ones([1,size(ChargeLabelsInternalUniqueBottom{1},2)]),+ones([1,size(ChargeLabelsInternalUniqueTop{1},2)])];
            
            %{
            if SymHandle.IsNonAbelian
                
                %{
                SideStructure = zeros(size(Structure));
                [~,IndexExt] = sort(-Structure(Structure<0),'ascend');
                [~,IndexInt] = sort(Structure(Structure>0),'ascend');
                
                SideStructure(Structure<0) = ChargeSide(IndexExt);
                SideStructure(Structure>0) = -ChargeSideInt(IndexInt);
                
                PowerTerms = zeros(size(Structure));
                OutSameSide = SideStructure(1,:)==SideStructure(2,:);
                PowerTerms = 1*repmat(OutSameSide,[2,1]);
                PowerTerms(:,2:end) = PowerTerms(:,2:end) + (~OutSameSide([1,1],2:end)).*(1-2*([ChargeSideInt;ChargeSideInt] ~= SideStructure(:,2:end) ));
                
                PowerTermsInt = PowerTerms(Structure>0); PowerTermsInt = PowerTermsInt(:)';
                PowerTermsInt = PowerTermsInt(IndexInt) + (-1).^OutSameSide(:,2:end);
                
                PowerTermsExt = PowerTerms(Structure<0); PowerTermsExt = PowerTermsExt(:)';
                PowerTermsExt = PowerTermsExt(IndexExt);
                
                Dim = SymHandle.Dim;
                DimCorrections = prod(reshape(Dim([TempChargeLabelsExternal,TempChargeLabelsInternal],2),size([TempChargeLabelsExternal,TempChargeLabelsInternal]))...
                    .^(repmat([PowerTermsExt,PowerTermsInt],[size(TempChargeLabelsExternal,1),1])/4),2);
                %}
                DimCorrections = prod(reshape(Dim(TempChargeLabelsExternal,2),size(TempChargeLabelsExternal)),2).^(-1/4);
                
                for kk = 1:numel(TempTensorEntries)
                    TempTensorEntries{kk} = TempTensorEntries{kk}*DimCorrections(kk);
                end
            end
            %}
            
            
            Tensor = SymTensor(TempTensorEntries, Structure, ChargeDirections, TempMultiplicitiesInternal, TempChargeLabelsInternal,...
                TempChargeLabelsExternal, Braidings, BraidingDirections,[],ChargeLegDimensions,SymHandle,ChargeSide, ChargeSideInt);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Restructure if we wanted to for some reason
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isempty(ReshapeData) && ~isequal(ReshapeData,0)
                
                %don't need to worry about this yet
                
%                Tensor = Tensor.SymReshape(InverseOfReshapeData);
                
            end
            
            
        end
        
        function RotationModifications = FromStandard(Tensors, RotateTensor, Invert, Adjust)
            
            if isa(Tensors, 'SymTensor')
                Tensors = {Tensors};
                FlagSingleOutput = true;
            else
                if iscell(Tensors)
                    if ~isa(Tensors{1},'SymTensor')&& ~iscell(Tensors{1})
                        FlagSingleOutput = true;
                        Tensors = {Tensors};
                    else
                        FlagSingleOutput = false;
                    end
                else
                    FlagSingleOutput = false;
                end
            end
            
            if ~iscell(Tensors)
                error('Error: Tensors Should be a cell')
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            FlagTensorInput = false([numel(Tensors),1]);
            for ll = 1:numel(Tensors)
                if isa(Tensors{ll}, 'SymTensor')
                    FlagTensorInput(ll) = true;
                else
                    FlagTensorInput(ll) = false;
                    %do a check that we have everything
                end
            end
            
            
            if nargin<2
                RotateTensor = ones([1,numel(Tensors)]);
            elseif isempty(RotateTensor)
                RotateTensor = ones([1,numel(Tensors)]);
            end
            
            if nargin<3
                Invert = false;
            end
            
            if nargin<4
                Adjust = true;
            end
            
            RotationModifications = cell(size(Tensors));
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now All the prep work has been done
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            
            
            
            
            
            
            
            for kk = 1:numel(Tensors)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Construct ChargeSide, Structure, ChargeDirections and
                %CyclicNumbers from our details
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if FlagTensorInput(kk)
                    ChargeSide = Tensors{kk}.ChargeSide;
                    StructureFull = Tensors{kk}.Structure;
                    ChargeDirectionsFull = Tensors{kk}.ChargeDirections;
                    ChargeSideInt = Tensors{kk}.ChargeSideInt;
                    CyclicNumbers = Tensors{kk}.CyclicLegsBraided;
                    
                else
                    StructureFull = Tensors{kk}{1};
                    ChargeDirectionsFull = Tensors{kk}{2};
                    ChargeSide = Tensors{kk}{3};
                    ChargeSideInt = Tensors{kk}{4};
                    
                    %if We don't have CyclicNumbers here assume that we are
                    %counting from left to right on top and bottom sides
                    if numel(Tensors{kk})<5
                        Numbers = 1:numel(ChargeSide);
                        CyclicNumbers = Numbers(ChargeSide == -1);
                        CyclicNumbers = [Numbers(ChargeSide == +1), CyclicNumbers(end:-1:1)];
                    else
                        CyclicNumbers = Tensors{kk}{5};
                    end
                    
                end

                RotationModifications{kk} = cell([0,2]);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Perform Edge cases first
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if isempty(StructureFull)
                    
                    %case of either 0 or 1 legs
                    if isempty(ChargeSide)
                        %no legs
                        %do nothing
                        continue;
                        
                    else
                        %one leg
                        
                        if ChargeSide == -1;
                            RotationModifications{kk} = {'RCCWInv',[-1,nan,0,1;ChargeDirections,0,+1,0]};
                        end
                        
                        continue;
                        
                    end
                    
                    
                elseif size(StructureFull,2) == 1
                    
                    %then we have 2 legs
                    
                    if ChargeSide(1) == +1 && ChargeSide(2) == +1
                        aLabel = -CyclicNumbers(1);
                        aDir = ChargeDirectionsFull(StructureFull == aLabel);
                        bLabel = -CyclicNumbers(2);
                        bDir = ChargeDirectionsFull(StructureFull == bLabel);
                    elseif ChargeSide(1) == -1 && ChargeSide(2) == -1
                        aLabel = -CyclicNumbers(2);
                        aDir = ChargeDirectionsFull(StructureFull == aLabel);
                        bLabel = -CyclicNumbers(1);
                        bDir = ChargeDirectionsFull(StructureFull == bLabel);
                    else
                        aLabel = -find(ChargeSide == +1);
                        aDir = ChargeDirectionsFull(StructureFull == aLabel);
                        bLabel = -find(ChargeSide == -1);
                        bDir = ChargeDirectionsFull(StructureFull == bLabel);
                        if ~isequal(CyclicNumbers,[-aLabel,-bLabel])
                            error('Error: CyclicLegs needs to be clockwise from left top, for this two leg Tensor it is not of that form')
                        end
                    end
                    
                    if ChargeSide(1) == +1 && ChargeSide(2) == +1
                        if RotateTensor(kk) == 1
                            %do nothing
                            First = CyclicNumbers(1);
                        else
                            RotationModifications{kk} = {'pDown',[aLabel;aDir]};
                            First = CyclicNumbers(2);
                        end
                        
                        %renaming if needed
                        if First~=1
                            RotationModifications{kk} = [RotationModifications{kk};{'Swap',[-2,-1;0,0;0,0;-2,-1;0,0;0,0]}];
                        end
                    elseif ChargeSide(1) == -1 && ChargeSide(2) == -1
                        if RotateTensor(kk)==1 
                            RotationModifications{kk} = {'pDown',[aLabel;-aDir]};
                            First = CyclicNumbers(1);
                        elseif ~Adjust
                            %then do nothing
                            First = CyclicNumbers(2);
                        else
                            RotationModifications{kk} = {'pDown',[aLabel;-aDir]};
                            First = CyclicNumbers(2);
                        end
                        
                        %renaming if needed
                        if First~=1
                            RotationModifications{kk} = [RotationModifications{kk};{'Swap',[-2,-1;0,0;0,0;-2,-1;0,0;0,0]}];
                        end
                    else
                        if RotateTensor(kk) == 1
                            if Adjust
                                RotationModifications{kk} = [{'pDown';'RCWInv'},...
                                    {[bLabel;-bDir];[aLabel, bLabel, NaN, NaN;aDir,bDir,0,0]}];
                                First = -aLabel;
                            else
                                RotationModifications{kk} = [{'RCWInv'},...
                                    {[aLabel, bLabel, NaN, NaN;aDir,bDir,0,0]}];
                                First = -aLabel;
                            end
                        else
                            RotationModifications{kk} = [{'RCCWInv'},...
                                {[bLabel, aLabel, NaN, NaN;bDir,aDir,0,0]}];
                            First = -bLabel;
                        end
                        
                        %renaming if needed
                        if First~=1
                            RotationModifications{kk} = [RotationModifications{kk};{'Swap',[-2,-1;0,0;0,0;-2,-1;0,0;0,0]}];
                        end
                    end
                    
                else %if size(StructureFull,2)>1
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %zeroth we want to make sure the StructureFull is of the
                %form that bottom one indicates most right or most bottom (in
                %that order) Using ChargeSide StructureFull and
                %CyclicNumbers
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                FixedStructure = false([1,size(StructureFull,2)]);
                
                ChargeSidetructure = zeros(size(StructureFull));
                
                for pp = 1:(size(StructureFull,2)+1)
                    ChargeSidetructure(StructureFull == -pp) = ChargeSide(pp);
                    ChargeDirectionExt(pp) = ChargeDirectionsFull(StructureFull==-pp);
                end
                for pp = 1:(size(StructureFull,2)-1)
                    ChargeSidetructure(StructureFull == pp) = ChargeSideInt(pp);
                    ChargeDirectionInt(pp) = ChargeDirectionsFull(StructureFull==pp);
                end
                
                if any(any(ChargeSidetructure == 0))
                    error('Affirmation Error: I should have filled structure in completely')
                end
                
                ChargeSidetructureSame = ChargeSidetructure(1,:) == ChargeSidetructure(2,:);
                
                ModChargeSide = ChargeSide(CyclicNumbers);
                TopChargesOriginal = CyclicNumbers(ModChargeSide == +1);
                BottomChargesOriginal = CyclicNumbers(ModChargeSide == -1);
                
                TopChargesModified = -TopChargesOriginal;
                BottomChargesModified = -BottomChargesOriginal(end:-1:1);
                
                ChargeSidetructureAccessible = StructureFull<0;
                
                WorkWith = find(ChargeSidetructureAccessible(1,:)&ChargeSidetructureAccessible(2,:));
                
                while ~isempty(WorkWith)
                    
                    for pp = WorkWith
                        ChargeSidetructureAccessible(StructureFull == pp-1) = true;
                        
                        if ChargeSidetructureSame(pp)
                            if ChargeSidetructure(1,pp) == +1
                                FirstLoc = find(TopChargesModified == StructureFull(1,pp));
                                SecondLoc = find(TopChargesModified == StructureFull(2,pp));
                                
                                if FirstLoc<SecondLoc
                                    %then things are safe
                                    
                                    if FirstLoc ~= SecondLoc-1;
                                        error('Error: This structure doesn''t obey the given cyclic order')
                                    end
                                else
                                    %then things aren't safe
                                    StructureFull([1,2],pp) = StructureFull([2,1],pp);
                                    ChargeDirectionsFull([1,2],pp) = ChargeDirectionsFull([2,1],pp);
                                    
                                    if FirstLoc ~= SecondLoc+1;
                                        error('Error: This structure doesn''t obey the given cyclic order')
                                    end
                                    
                                end
                                
                                if pp ~=1
                                    if ChargeSideInt(pp-1)~=+1
                                        error('Error: ChargeSideInt and StructureFull give a vertex which has nothing coming out the bottom in structure')
                                    end
                                end
                                
                                TopChargesModified(FirstLoc) = pp-1;
                                TopChargesModified(SecondLoc) = [];
                                FixedStructure(pp) = true;
                                
                            else
                                FirstLoc = find(BottomChargesModified == StructureFull(1,pp));
                                SecondLoc = find(BottomChargesModified == StructureFull(2,pp));
                                
                                if FirstLoc<SecondLoc
                                    %then things are safe
                                    
                                    if FirstLoc ~= SecondLoc-1;
                                        error('Error: This structure doesn''t obey the given cyclic order')
                                    end
                                else
                                    %then things aren't safe
                                    StructureFull([1,2],pp) = StructureFull([2,1],pp);
                                    ChargeDirectionsFull([1,2],pp) = ChargeDirectionsFull([2,1],pp);
                                    
                                    if FirstLoc ~= SecondLoc+1;
                                        error('Error: This structure doesn''t obey the given cyclic order')
                                    end
                                    
                                end
                                
                                if pp ~=1
                                    if ChargeSideInt(pp-1)~=-1
                                       error('Error: ChargeSideInt and StructureFull give a vertex which has nothing coming out the top in structure')
                                    end
                                end
                                
                                BottomChargesModified(FirstLoc) = pp-1;
                                BottomChargesModified(SecondLoc) = [];
                                FixedStructure(pp) = true;
                                
                            end
                        else
                            if any(TopChargesModified == StructureFull(1,pp))
                                if any(BottomChargesModified == StructureFull(2,pp))
                                    %all good
                                else
                                    %not possible
                                    error('Error: StructureFull does not agree with ChargeSide/ChargeSideInt');
                                end
                            elseif any(BottomChargesModified == StructureFull(1,pp))
                                if any(TopChargesModified == StructureFull(2,pp))
                                    %need to swap
                                    StructureFull([1,2],pp) = StructureFull([2,1],pp);
                                    ChargeDirectionsFull([1,2],pp) = ChargeDirectionsFull([2,1],pp);
                                else
                                    %not possible
                                    error('Error: StructureFull does not agree with ChargeSide/ChargeSideInt');
                                end
                            else
                                error('Error: StructureFull does not agree with ChargeSide/ChargeSideInt');
                            end
                            
                            
                            FirstLoc = find(TopChargesModified == StructureFull(1,pp));
                            SecondLoc = find(BottomChargesModified == StructureFull(2,pp));
                            if pp ~= 1
                                
                            if FirstLoc == 1 && SecondLoc == 1
                                if ChargeSideInt(pp-1) == +1;
                                    TopChargesModified(1) = pp-1;
                                    BottomChargesModified(1) = [];
                                else
                                    TopChargesModified(1) = [];
                                    BottomChargesModified(1) = pp-1;
                                end
                            elseif FirstLoc == numel(TopChargesModified) && SecondLoc == numel(BottomChargesModified)
                                if ChargeSideInt(pp-1) == +1;
                                    TopChargesModified(end) = pp-1;
                                    BottomChargesModified(end) = [];
                                else
                                    TopChargesModified(end) = [];
                                    BottomChargesModified(end) = pp-1;
                                end
                            else
                                error('Error: StructureFull does not agree the given cyclic order');
                            end
                            else
                                if numel(TopChargesModified)~=1 || numel(BottomChargesModified)~=1
                                    error('Error: StructureFull does not work');
                                end
                            end
                            
                            FixedStructure(pp) = true;
                            
                        end
                        
                    end
                    
                    WorkWith = find(ChargeSidetructureAccessible(1,:)&ChargeSidetructureAccessible(2,:)&~FixedStructure);
                end
                
                if any(~FixedStructure)
                    error('Error: StructureFull does not connect in the correct manner')
                end
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %first we work out the permutation and rotation effects:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                ChargeSideCyclic = ChargeSide(CyclicNumbers);
                TopNumbers = CyclicNumbers(ChargeSideCyclic == +1);
                BottomNumbers = CyclicNumbers(ChargeSideCyclic == -1); %these are already right to left
                %BottomNumbers = BottomNumbers(end:-1:1);
                LegType = zeros([1,numel(CyclicNumbers)]);
                
                if RotateTensor(kk) == 1
                    LegType(TopNumbers) = 1; % Good
                    LegType(BottomNumbers) = 4; %UpRight
                    
                    if Adjust
                        for ll = BottomNumbers
                            ChargeDir = ChargeDirectionsFull(StructureFull == -ll);
                            if numel(ChargeDir) ~= 1
                                error('Affirmation Error: there are too many charge directions for this pDown');
                            end
                            RotationModifications{kk} = [RotationModifications{kk};{'pDown',[-ll;-ChargeDir]}];
                        end
                    end
                    
                    LegsOutOrder = [TopNumbers, BottomNumbers];
                elseif RotateTensor(kk)>numel(TopNumbers)
                    %then we don't do any rotations on the top
                    LegType(TopNumbers) = 1; %Good
                    LegType(BottomNumbers(1:(RotateTensor(kk)-numel(TopNumbers)-1))) = 4; %UpRight
                    if Adjust
                        for ll = BottomNumbers(1:(RotateTensor(kk)-numel(TopNumbers)-1))
                            ChargeDir = ChargeDirectionsFull(StructureFull == -ll);
                            if numel(ChargeDir) ~= 1
                                error('Affirmation Error: there are too many charge directions for this pDown');
                            end
                            RotationModifications{kk} = [RotationModifications{kk};{'pDown',[-ll;-ChargeDir]}];
                        end
                    end
                    LegType(BottomNumbers((RotateTensor(kk)-numel(TopNumbers)):end)) = 3; %UpLeft
                    
                    LegsOutOrder = [BottomNumbers((RotateTensor(kk)-numel(TopNumbers)):end), TopNumbers, BottomNumbers(1:(RotateTensor(kk)-numel(TopNumbers)-1))];
                else
                    %then we have all bottoms are upLeft
                    LegType(BottomNumbers) = 3; %UpLeft
                    LegType(TopNumbers(RotateTensor(kk):end)) = 2; %Rotate
                    
                    %for ll = TopNumbers(RotateTensor(kk):end)
                    %    ChargeDir = ChargeDirectionsFull(StructureFull == -ll);
                    %    if numel(ChargeDir) ~= 1
                    %        error('Affirmation Error: there are too many charge directions for this pDown');
                    %    end
                    %    RotationModifications{kk} = [RotationModifications{kk};{'Round',[-ll;ChargeDir]}];
                    %end
                    
                    LegType(TopNumbers(1:(RotateTensor(kk)-1))) = 1; %Good
                    
                    LegsOutOrder = [TopNumbers(RotateTensor(kk):end), BottomNumbers, TopNumbers(1:(RotateTensor(kk)-1))];
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %second thing we want to do is to update the verticies
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                VertexType = zeros([1,size(StructureFull,2)]);
                VertexUsed = false(size(VertexType));
                InternalLegsUseable = false([1,numel(CyclicNumbers)-2]);
                
                WorkingTop = -TopNumbers;
                WorkingBottom = -BottomNumbers(end:-1:1); %this is left to right
                StructureTemp = StructureFull;
                StructureUsable = StructureFull<0;
                ChargeDirectionsTemp = ChargeDirectionsFull;
                
                InternalInverse = zeros(size(VertexType));
                InternalSpecial = false(size(VertexType));
                
                while any(~VertexUsed)
                    
                    LocationsUpdatable = find((~VertexUsed)&all(StructureUsable,1));
                    
                    if isempty(LocationsUpdatable)
                        error('Affirmation Error: We can''t finish updating the verticies because apparently there are no more valid ones to update while there are unupdated verticies')
                    elseif any(LocationsUpdatable == 1)
                        if numel(LocationsUpdatable)~=1
                            error('Affirmation Error: the core can''t be updatable at the same time as any other vertex')
                            %if we get the core vertex then we should have
                            %updated everything else already.
                        end
                    end
                    
                    for ll = LocationsUpdatable
                        
                        
                        NumberOfInverses = 0;
                        SpecialInverse = false;
                        InverseSide = 0;
                        
                        aLabel = StructureTemp(1,ll);
                        if aLabel<0
                            aType = LegType(-aLabel);
                        else
                            aType = LegType(numel(CyclicNumbers)+aLabel+1);
                            NumberOfInverses = NumberOfInverses+InternalInverse(aLabel);
                            %SpecialInverse = SpecialInverse || InternalSpecial(aLabel);
                            %if InternalInverse(aLabel)
                            %    InverseSide = 1;
                            %end
                        end
                        aDirection = ChargeDirectionsTemp(1,ll);
                        
                        
                        bLabel = StructureTemp(2,ll);
                        if bLabel<0
                            bType = LegType(-bLabel);
                        else
                            bType = LegType(numel(CyclicNumbers)+bLabel+1);
                            NumberOfInverses = NumberOfInverses+InternalInverse(bLabel);
                            %SpecialInverse = SpecialInverse || InternalSpecial(bLabel);
                            %if InternalInverse(bLabel)
                            %    InverseSide = 2;
                            %end
                        end
                        bDirection = ChargeDirectionsTemp(2,ll);
                        
                        if ll ~= 1
                            cDirection = ChargeDirectionsTemp(StructureTemp==(ll-1));
                            cSide = ChargeSideInt(ll-1);
                            cLabel = ll-1;
                            MultLabel = ll-1;
                        else
                            cLabel = NaN;
                            MultLabel = NaN;
                            cDirection = 0;
                            cSide = -1 + 2*(aType<3 || bType<3);
                        end
                        
                        
                        
                        
                        %Now work out the cases:
                        
                        VertexUsed(ll) = true;
                        if aType<3 && bType<3 % both are top type
                            
                            
                            
                            
                            
                            %now actually work out everything
                            
                            if aType == 1 && bType == 1 %both good
                                LegType(numel(CyclicNumbers)+ll) = 1; %good
                                VertexType(ll) = 1; %good
                                CellUpdate = cell([0,2]);
                                %disp('1,1')
                            elseif aType == 2 && bType == 2 %both rotated
                                VertexType(ll) = 1; %good
                                LegType(numel(CyclicNumbers)+ll) = 2; %rotated
                                %disp('2,2')
                                CellUpdate = [{'RCW';'RCCWInv';'RCW';'RCCWInv';'RCW';'RCCWInv'},...
                                    {[aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0];...
                                    [cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                    [cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                    [bLabel,cLabel,aLabel,MultLabel,MultLabel;bDirection,-cDirection,-aDirection,0,0];...
                                    [bLabel,cLabel,aLabel,MultLabel,MultLabel;bDirection,-cDirection,-aDirection,0,0];...
                                    [aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0]}];
                                %RotationModifications{kk} = [RotationModifications{kk};CellUpdate];
                            else
                                if ~(aType == 1 && bType == 2)
                                    error('Affirmation Error: A good/Rotated vertex doesn''t agree with what I''m expecting');
                                end
                                VertexType(ll) = 2; %inverted
                                LegType(numel(CyclicNumbers)+ll) = 2; %rotated
                                %disp('1,2')
                                CellUpdate = [{'RCW';'RCCWInv';'RCW'},{[aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0];...
                                    [cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                    [cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0]}];
                                %RotationModifications{kk} = [RotationModifications{kk};CellUpdate];
                                
                                
                                %Swap a and b
                                ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                            end
                            
                        elseif aType>2 && bType>2
                            %both are bottom type
                            
                            
                            
                            %now actually work out everything
                            
                            if aType == 3 && bType == 3 %both UpLeft
                                VertexType(ll) = 1; %good
                                LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                %disp('3,3')
                                CellUpdate = [{'RCCWInv';'RCW';'RCCWInv';'pDownInv'},{[aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0];...
                                    [aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0];...
                                    [bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0];[cLabel;+cDirection]}];
                                
                                
                                
                                %Swap a and b
                                ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                                 
                            elseif aType == 4 && bType == 4 %both UpRight
                                VertexType(ll) = 1; %good
                                LegType(numel(CyclicNumbers)+ll) = 4; %UpRight
                                %disp('4,4')
                                CellUpdate = [{'RCWInv';'pDown';'RCCW';'RCWInv'},{[cLabel,bLabel,aLabel,MultLabel,MultLabel;-cDirection,bDirection,-aDirection,0,0];...
                                    [cLabel;-cDirection];[cLabel,bLabel,aLabel,MultLabel,MultLabel;-cDirection,bDirection,-aDirection,0,0];...
                                    [bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0]}];
                                
                                
                                
                                %Swap a and b
                                ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                                
                            else
                                if ~(aType == 3 && bType == 4)
                                    error('Affirmation Error: An UpLeft/UpRight vertex doesn''t agree with what I''m expecting');
                                end
                                VertexType(ll) = 2; %inverted
                                LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                %disp('3,4')
                                CellUpdate = [{'pDownInv'},{[cLabel;+cDirection]}];
                                
                            end
                            
                            
                            
                        else
                            
                            %now we know that there is one top and one bottom
                            %leg:
                            
                            if aType < 3 % if a on Top
                                %Do nothing, everything is fine
                                
                            else
                                %in this case we should swap
                                error('AffirmationError: it appears the flip did not happen problem for tensor %i, at vertex %i', kk,ll)
                                
                            end
                            
                            
                            if ll~=1
                                FlagCoreToRight = any(StructureFull(1,:) == ll-1);
                            else
                                FlagCoreToRight = true;
                            end 
                            
                            if cSide == -1 
                                FlagCoreComingFromAbove = true;
                            else
                                FlagCoreComingFromAbove = false;
                            end
                            
                            
                            %now do all the checks as above:
                            if ~FlagCoreToRight
                                if ~FlagCoreComingFromAbove
                                    %now we know that this is far right (so core is to the left) and the core comes into the vertex from below 
                                    
                                    
                                    if aType == 1 && bType == 4 %UpRight and Good
                                        VertexType(ll) = 1; %good
                                        LegType(numel(CyclicNumbers)+ll) = 1; %good
                                        disp('CBL 1,4')
                                        CellUpdate = [{'RCWInv'},{[aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0]}];
                                        
                                    elseif aType == 1 && bType == 3 %UpLeft and Good
                                        VertexType(ll) = 2; %inverse
                                        LegType(numel(CyclicNumbers)+ll) = 2; %rotated
                                        disp('CBL 1,3')
                                        CellUpdate = [{'RCCWInv';'RCW'},{[cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                            [cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0]}];
                                        
                                        %Swap a and b
                                        ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                        StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                        StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                                    else
                                        if ~(aType == 2 && bType == 3)
                                            error('Affirmation Error: An unallowed combination appears on a vertex');
                                        end
                                        VertexType(ll) = 1; %good
                                        LegType(numel(CyclicNumbers)+ll) = 2; %rotated
                                        disp('CBL 2,3')
                                        CellUpdate = [{'RCCWInv';'RCW';'RCCWInv';'RCW';'RCCWInv'},...
                                            {[cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                            [cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                            [bLabel,cLabel,aLabel,MultLabel,MultLabel;bDirection,-cDirection,-aDirection,0,0];...
                                            [bLabel,cLabel,aLabel,MultLabel,MultLabel;bDirection,-cDirection,-aDirection,0,0];...
                                            [aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0]}];
                                        
                                    end
                                    
                                    
                                    
                                    
                                    
                                    
                                else
                                    %now we know that this is far right (so core is to the left) and the core comes into the vertex from above 
                                    
                                    
                                    if aType == 1 && bType == 4 %UpRight and Good
                                        VertexType(ll) = 1; %good
                                        LegType(numel(CyclicNumbers)+ll) = 4; %UpRight
                                        disp('CAL 1,4')
                                        CellUpdate = [{'RCCW';'RCWInv';'pDown'},{[cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                            [aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0];[cLabel;-cDirection]}];
                                        
                                    elseif aType == 1 && bType == 3 %UpLeft and Good
                                        VertexType(ll) = 2; %inverse
                                        LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                        disp('CAL 1,3')
                                        CellUpdate = [{'RCWInv';'pDown'},{[cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                            [cLabel;+cDirection]}];
                                        
                                        %Swap a and b
                                        ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                        StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                        StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                                        
                                    else
                                        if ~(aType == 2 && bType == 3)
                                            error('Affirmation Error: An unallowed combination appears on a vertex');
                                        end
                                        VertexType(ll) = 1; %good
                                        LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                        disp('CAL 2,3')
                                        CellUpdate = [{'RCW';'RCCWInv';'RCW';'RCCWInv';'pDownInv'},{[cLabel,aLabel,bLabel,MultLabel,MultLabel;-cDirection,aDirection,-bDirection,0,0];...
                                            [bLabel,cLabel,aLabel,MultLabel,MultLabel;bDirection,-cDirection,-aDirection,0,0];...
                                            [bLabel,cLabel,aLabel,MultLabel,MultLabel;bDirection,-cDirection,-aDirection,0,0];...
                                            [aLabel,bLabel,cLabel,MultLabel,MultLabel;aDirection,bDirection,cDirection,0,0];[cLabel;+cDirection]}];
                                        
                                    end
                                    
                                    
                                    
                                end
                            else
                                if ~FlagCoreComingFromAbove
                                    %now we know that this is far left (so core is to the right) and the core comes into the vertex from below 
                                    
                                    
                                    
                                    if aType == 1 && bType == 4 %UpRight and Good
                                        VertexType(ll) = 2; %inverted
                                        LegType(numel(CyclicNumbers)+ll) = 1; %good
                                        UseSpecial = true;
                                        %disp('CBR 1,4')
                                        CellUpdate = [{'RCWInv';'RCCW'},{[aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0];...
                                            [aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0]}];
                                        
                                    elseif aType == 1 && bType == 3 %UpLeft and Good
                                        VertexType(ll) = 1; %good
                                        LegType(numel(CyclicNumbers)+ll) = 1; %good
                                        %disp('CBR 1,3')
                                        CellUpdate = [{'RCCWInv'},{[bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0]}];
                                        
                                        %Swap a and b
                                        ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                        StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                        StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                                        
                                    else
                                        if ~(aType == 2 && bType == 3)
                                            error('Affirmation Error: An unallowed combination appears on a vertex');
                                        end
                                        VertexType(ll) = 2; %inverse
                                        LegType(numel(CyclicNumbers)+ll) = 2; %rotated
                                        %disp('CBR 2,3')
                                        CellUpdate = [{'RCCWInv';'RCW';'RCCWInv';'RCW'},...
                                            {[bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0];...
                                            [bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0];...
                                            [cLabel,bLabel,aLabel,MultLabel,MultLabel;-cDirection,bDirection,-aDirection,0,0];...
                                            [cLabel,bLabel,aLabel,MultLabel,MultLabel;-cDirection,bDirection,-aDirection,0,0]}];
                                        
                                    end
                                    
                                    
                                else
                                    %now we know that this is far left (so core is to the right) and the core comes into the vertex from above 
                                    
                                    
                                    if aType == 1 && bType == 4 %UpRight and Good
                                        VertexType(ll) = 2; %inverted
                                        LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                        disp('CAR 1,4')
                                        CellUpdate = [{'RCCW';'pDownInv'},{[aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0];...
                                            [cLabel;+cDirection]}];
                                        
                                        
                                    elseif aType == 1 && bType == 3 %UpLeft and Good
                                        VertexType(ll) = 1; %good
                                        LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                        disp('CAR 1,3')
                                        CellUpdate = [{'RCW';'RCCWInv';'pDownInv'},{[aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0];...
                                            [bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0];[cLabel;+cDirection]}];
                                        
                                        %Swap a and b
                                        ChargeDirectionsTemp([1,2],ll) = ChargeDirectionsTemp([2,1],ll);
                                        StructureTemp([1,2],ll) = StructureTemp([2,1],ll);
                                        StructureUsable([1,2],ll) = StructureUsable([2,1],ll);
                                        
                                    else
                                        if ~(aType == 2 && bType == 3)
                                            error('Affirmation Error: An unallowed combination appears on a vertex');
                                        end
                                        VertexType(ll) = 2; %inverted
                                        LegType(numel(CyclicNumbers)+ll) = 3; %UpLeft
                                        UseSpecial = true;
                                        
                                        error('Error: CodeHERA')
                                        
                                        CellUpdate = [{'RCW';'RCCWInv';'RCW';'RCCWInv';'RCW'},{[aLabel,cLabel,bLabel,MultLabel,MultLabel;aDirection,-cDirection,-bDirection,0,0];...
                                            [bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0];...
                                            [bLabel,aLabel,cLabel,MultLabel,MultLabel;bDirection,aDirection,cDirection,0,0];...
                                            [cLabel,bLabel,aLabel,MultLabel,MultLabel;-cDirection,bDirection,-aDirection,0,0];...
                                            [cLabel,bLabel,aLabel,MultLabel,MultLabel;-cDirection,bDirection,-aDirection,0,0]}];
                                        
                                    end
                                    
                                end
                            end
                            
                        end
                        
                        StructureUsable(StructureTemp == (ll-1)) = true;
                        RotationModifications{kk} = [RotationModifications{kk};CellUpdate];
                        
                    end
                    
                end
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %at this point we can make the assumption that all good
                %verticies and inverted verticies have labelings a on left and
                %b on right.
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %third step if to convert all inverse verticies into
                %good verticies
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                
                if VertexType(1) == 2
                    %in this case then the core is already a inverse
                    %vertex. This we can immediently correct:
                    
                    %aLabel and bLabel are set as the last thing we
                    %need to check is the core:
                    
                    %CellUpdate = [{'RCWInv';'RCCWInv'}, {[NaN, bLabel, aLabel, NaN, NaN; 0, bDirection, aDirection, 0,0];[aLabel, bLabel, NaN, NaN, NaN; +aDirection, +bDirection, 0,0,0]}];
                    
                    %RotationModifications{kk} = [RotationModifications{kk};CellUpdate];
                    
                    %CellUpdate = [{'pDownInv'}, {[ aLabel; aDirection]}];
                    VertexType(1) = 1;%it is now good.
                end
                %RotationModifications{kk} = [RotationModifications{kk};CellUpdate];
                Loc = find(VertexType == 2)-1;
                
                if ~isempty(Loc)
                    if any(Loc<1)
                        error('Affirmation Error: We shouldn''t get that the core is a goal in part 3')
                    end
                end
                
                SwapDataOriginal = [-(1:numel(ChargeSide)),1:numel(ChargeSideInt);zeros([1,numel(ChargeSide)+numel(ChargeSideInt)]);1:numel(ChargeSideInt),zeros([1,numel(ChargeSide)])];
                SwapDataOriginal = [SwapDataOriginal;SwapDataOriginal];
                
                Counter = 1;
                FlagBegining = false;
                while ~isempty(Loc)
                    if Counter>numel(Loc)
                        if FlagBegining == true
                            error('Affirmation Error: I will be getting an infinite loop in part 3')
                        end
                        Counter = 1;
                        FlagBeginning = true;
                    end
                    
                    
                    %set the Location I'm using, note that we aren't changing
                    %anything else in Loc as we don't allow ourselves to move
                    %through inverses
                    %create a path which starts at our goal.
                    %create an empty list which stores if the connection is the
                    %first or second entry (this is smaller then path as we
                    %don't need anything at our goal
                    LocUsing = Loc(Counter);
                    Path = LocUsing;
                    IsTopPath = false([0,1]);
                    
                    while LocUsing ~=0
                        LocUsing = find(any(StructureTemp == LocUsing,1))-1;
                        Path = [Path, LocUsing];
                        if isempty(LocUsing)
                            error('Affirmation Error: We couldn''t find the previous vertex from this Loc Using')
                        elseif numel(LocUsing)>1
                            error('Affirmation Error: We found too many previous verticies, we should only have found one')
                        end
                        if StructureTemp(1,LocUsing+1) == Path(end-1)
                            IsTopPath = [IsTopPath,true];
                        else
                            IsTopPath = [IsTopPath,false];
                        end
                    end
                    
                    if sum(VertexType(Path+1) == 2) >1
                        Counter = Counter+1;
                        continue;
                    end
                    
                    
                    
                    %now that we have worked out the path we start with the
                    %core and move through to the point before the goal.
                    
                    CounterPath = 0;
                    AdjustPath = cell([0,2]);
                    
                    for ll = (length(Path)-1):-1:2
                        %we don't want to do the first entry as we aren't
                        %moving the core to there
                        CounterPath = CounterPath+1;
                        %Note that Path(end) = 0
                        
                        if IsTopPath(ll)
                            LocL = [1,2*Path(ll)+1, 2*Path(ll)+2,2];
                            LocL2 = LocL([2,3,4,1]);
                            
                            dLabel = StructureTemp(LocL(1));
                            aLabel = StructureTemp(LocL(2));
                            bLabel = StructureTemp(LocL(3));
                            cLabel = StructureTemp(LocL(4));
                            
                            dDir = ChargeDirectionsTemp(LocL(1));
                            aDir = ChargeDirectionsTemp(LocL(2));
                            bDir = ChargeDirectionsTemp(LocL(3));
                            cDir = ChargeDirectionsTemp(LocL(4));
                            
                            StructureTemp(LocL) = StructureTemp(LocL2);
                            SwapForward = ChargeDirectionsTemp(LocL(2)) == +1;
                            SwapBackward = ChargeDirectionsTemp(LocL(4)) == +1;
                            
                            ChargeDirectionsTemp(LocL(1)) = +1;
                            ChargeDirectionsTemp(LocL) = ChargeDirectionsTemp(LocL2);
                            
                            SwapData = SwapDataOriginal;
                            SwapData(2,SwapData(1,:)==dLabel) = SwapForward;
                            SwapData(1,SwapData(1,:)==dLabel) = aLabel;
                            %need to add in multiplicities
                            
                            SwapData(5,SwapData(4,:)==dLabel) = SwapBackward;
                            SwapData(4,SwapData(4,:)==dLabel) = cLabel;
                            %need to add in multiplicities
                            
                            MultLabel = dLabel;
%                            AdjustPath = [AdjustPath;[{'Swap';'RCCW';'RCWInv';'pDown'},{SwapData;[dLabel, bLabel, cLabel, MultLabel,MultLabel; -dDir,bDir,-cDir,0,0];...
%                                [bLabel, cLabel, dLabel, MultLabel,MultLabel; bDir,cDir,dDir,0,0];[cLabel;-cDir]}]];
                            AdjustPath = [AdjustPath;[{'pDown';'RCCW';'RCWInv';'Swap'},{[cLabel;-cDir];[aLabel, bLabel, cLabel, MultLabel,MultLabel; aDir,bDir,-cDir,0,0];...
                                [bLabel, cLabel, aLabel, MultLabel,MultLabel; bDir,cDir,-aDir,0,0];SwapData}]];
                            %disp('AdjustedOther')
                            
                            
                        else
                            LocL = [1,2*Path(ll)+1, 2*Path(ll)+2,2];
                            LocL2 = LocL([4,1,2,3]);
                            
                            dLabel = StructureTemp(LocL(4));
                            aLabel = StructureTemp(LocL(1));
                            bLabel = StructureTemp(LocL(2));
                            cLabel = StructureTemp(LocL(3));
                            
                            dDir = ChargeDirectionsTemp(LocL(4));
                            aDir = ChargeDirectionsTemp(LocL(1));
                            bDir = ChargeDirectionsTemp(LocL(2));
                            cDir = ChargeDirectionsTemp(LocL(3));
                            
                            StructureTemp(LocL) = StructureTemp(LocL2);
                            SwapForward = ChargeDirectionsTemp(LocL(3)) == +1;
                            SwapBackward = ChargeDirectionsTemp(LocL(1)) == +1;
                            
                            ChargeDirectionsTemp(LocL(4)) = +1;
                            ChargeDirectionsTemp(LocL) = ChargeDirectionsTemp(LocL2);
                            
                            SwapData = SwapDataOriginal;
                            SwapData(2,SwapData(1,:)==dLabel) = SwapForward;
                            SwapData(1,SwapData(1,:)==dLabel) = cLabel;
                            %need to add in multiplicities
                            
                            SwapData(5,SwapData(4,:)==dLabel) = SwapBackward;
                            SwapData(4,SwapData(4,:)==dLabel) = aLabel;
                            %need to add in multiplicities
                            
                            MultLabel = dLabel;
                            
%                            AdjustPath = [AdjustPath;[{'pDownInv';'pDown';'RCW';'RCCWInv';'Swap'},{[cLabel;-cDir];[aLabel;aDir];...
%                                [bLabel, cLabel, dLabel, MultLabel,MultLabel; bDir,cDir,dDir,0,0];...
%                                [dLabel, bLabel, cLabel, MultLabel,MultLabel; -dDir,bDir,-cDir,0,0];SwapData}]];
%                            AdjustPath = [AdjustPath;[{'pDownInv';'RCW';'RCCWInv';'Swap'},{[cLabel;-cDir];...
%                                [bLabel, cLabel, dLabel, MultLabel,MultLabel; bDir,cDir,dDir,0,0];...
%                                [dLabel, bLabel, cLabel, MultLabel,MultLabel; -dDir,bDir,-cDir,0,0];SwapData}]];
                            AdjustPath = [AdjustPath;[{'pDownInv';'RCW';'RCCWInv';'Swap'},{[cLabel;-cDir];...
                                [bLabel, cLabel, aLabel, MultLabel,MultLabel; bDir,cDir,-aDir,0,0];...
                                [aLabel, bLabel, cLabel, MultLabel,MultLabel; aDir,bDir,-cDir,0,0];SwapData}]];
                            %disp('Adjusted')
                            
                        end
                    end
                    
                    VertexType(Path(1)+1) = 1;
                    Loc(Counter) = [];
                    FlagBeginning = false; %we have updated at least 1 Loc this time around.
                    
                    %special case: not shifting the core in the normal way
                    %here
                    
                    if ~IsTopPath(1)
                        error('Affirmation Error: the Inverse vertex should come out the left')
                    end
                    
                    LocL = [1,2*Path(1)+1, 2*Path(1)+2,2];
                    LocL2 = LocL([1,2,4,3]);
                    
                    dLabel = StructureTemp(LocL(4));
                    aLabel = StructureTemp(LocL(2));
                    bLabel = StructureTemp(LocL(3));
                    cLabel = StructureTemp(LocL(1));
                    
                    aDir = ChargeDirectionsTemp(LocL(2));
                    bDir = ChargeDirectionsTemp(LocL(3));
                    dDir = ChargeDirectionsTemp(LocL(4));
                    cDir = -bDir;
                    
                    StructureTemp(LocL) = StructureTemp(LocL2);
                    ChargeDirectionsTemp(LocL(1)) = -ChargeDirectionsTemp(LocL(3));
                    ChargeDirectionsTemp(LocL) = ChargeDirectionsTemp(LocL2);
                    
                    SwapData = SwapDataOriginal;
                    SwapData(2,SwapData(1,:)==cLabel) = 1;
                    SwapData(1,SwapData(1,:)==cLabel) = bLabel;
                    %need to add in multiplicities
                    
                    SwapData(5,SwapData(4,:)==cLabel) = 1;
                    SwapData(4,SwapData(4,:)==cLabel) = dLabel;
                    %need to add in multiplicities
                    
                    MultLabel = cLabel;
                    
                    AdjustPath = [AdjustPath;[{'Swap';'pDown';'pDownInv';'RCCWInv'},{SwapData;[dLabel;-dDir];[bLabel;-bDir];...
                        [aLabel, dLabel, cLabel, MultLabel,MultLabel; aDir,dDir,cDir,0,0]}]];
                    %disp('Random')
                    
                    RotationModifications{kk} = [RotationModifications{kk};AdjustPath];
                end
                
                %but further we want to make sure that the internal charges are
                %always going up (to set a convention and to make things easier
                %for part 4)
                
                FlipInternal = StructureTemp(StructureTemp>0 & ChargeDirectionsTemp == -1);
                ChargeDirectionsTemp(StructureTemp>0) = +1;
                
                if any(FlipInternal<0)
                    error('Affirmation Error: We are trying to change external legs')
                end
                
                if ~isempty(FlipInternal)
                    SwapData = SwapDataOriginal;
                    
                    for mm = FlipInternal'
                        SwapData(2,SwapData(1,:)==mm) = 1;
                        SwapData(5,SwapData(4,:)==mm) = 1;
                    end
                    RotationModifications{kk} = [RotationModifications{kk};{'Swap',SwapData}];
                end
                
                %From now one we can assume that the internal legs all point up
                %and all verticies are good
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %forth part is to get the final point on the right
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                LocUsing = find(any(StructureTemp==(-LegsOutOrder(end)),1))-1;
                if isempty(LocUsing)
                    error('Affirmation Error: we should be able to find at least one vertex which should be made right-most')
                end
                if numel(LocUsing)>1
                    error('Affirmation Error: there should only be one vertex which should be made right-most')
                end
                
                
                Path = LocUsing;
                IsTopPath = false([0,1]);
                
                while LocUsing ~=0
                    LocUsing = find(any(StructureTemp == LocUsing,1))-1;
                    Path = [Path, LocUsing];
                    if isempty(LocUsing)
                        error('Affirmation Error: We couldn''t find the previous vertex from this Loc Using')
                    elseif numel(LocUsing)>1
                        error('Affirmation Error: We found too many previous verticies, we should only have found one')
                    end
                    if StructureTemp(1,LocUsing+1) == Path(end)
                        IsTopPath = [IsTopPath,true];
                    else
                        IsTopPath = [IsTopPath,false];
                    end
                end
                
                
                
                
                %now that we have worked out the path we start with the
                %core and move through to the point before the goal.
                
                AdjustPath = cell([0,2]);
                for ll = (length(Path)-1):-1:1
                    %Note that Path(end) = 0
                    
                    if IsTopPath(ll)
                        LocL = [1,2*Path(ll)+1, 2*Path(ll)+2,2];
                        LocL2 = LocL([2,3,4,1]);
                        
                        dLabel = StructureTemp(LocL(1));
                        aLabel = StructureTemp(LocL(2));
                        bLabel = StructureTemp(LocL(3));
                        cLabel = StructureTemp(LocL(4));
                        
                        dDir = ChargeDirectionsTemp(LocL(1));
                        aDir = ChargeDirectionsTemp(LocL(2));
                        bDir = ChargeDirectionsTemp(LocL(3));
                        cDir = ChargeDirectionsTemp(LocL(4));
                        
                        StructureTemp(LocL) = StructureTemp(LocL2);
                        SwapForward = ChargeDirectionsTemp(LocL(2)) == +1;
                        SwapBackward = ChargeDirectionsTemp(LocL(4)) == +1;
                        
                        ChargeDirectionsTemp(LocL(1)) = +1;
                        ChargeDirectionsTemp(LocL) = ChargeDirectionsTemp(LocL2);
                        
                        SwapData = SwapDataOriginal;
                        SwapData(2,SwapData(1,:)==dLabel) = SwapForward;
                        SwapData(1,SwapData(1,:)==dLabel) = aLabel;
                        %need to add in multiplicities
                        
                        SwapData(5,SwapData(4,:)==dLabel) = SwapBackward;
                        SwapData(4,SwapData(4,:)==dLabel) = cLabel;
                        %need to add in multiplicities
                        
                        MultLabel = dLabel;
                        
                        AdjustPath = [AdjustPath;[{'pDown';'RCCW';'RCWInv';'Swap'},{[cLabel;-cDir];[aLabel, bLabel, cLabel, MultLabel,MultLabel; aDir,bDir,-cDir,0,0];...
                            [bLabel, cLabel, aLabel, MultLabel,MultLabel; bDir,cDir,-aDir,0,0];SwapData}]];
                        %disp('AdjustedOther')
                            
                    else
                        LocL = [1,2*Path(ll)+1, 2*Path(ll)+2,2];
                        LocL2 = LocL([4,1,2,3]);
                        
                        dLabel = StructureTemp(LocL(4));
                        aLabel = StructureTemp(LocL(1));
                        bLabel = StructureTemp(LocL(2));
                        cLabel = StructureTemp(LocL(3));
                        
                        dDir = ChargeDirectionsTemp(LocL(4));
                        aDir = ChargeDirectionsTemp(LocL(1));
                        bDir = ChargeDirectionsTemp(LocL(2));
                        cDir = ChargeDirectionsTemp(LocL(3));
                        
                        StructureTemp(LocL) = StructureTemp(LocL2);
                        SwapForward = ChargeDirectionsTemp(LocL(3)) == +1;
                        SwapBackward = ChargeDirectionsTemp(LocL(1)) == +1;
                        
                        ChargeDirectionsTemp(LocL(4)) = +1;
                        ChargeDirectionsTemp(LocL) = ChargeDirectionsTemp(LocL2);
                        
                        SwapData = SwapDataOriginal;
                        SwapData(2,SwapData(1,:)==dLabel) = SwapForward;
                        SwapData(1,SwapData(1,:)==dLabel) = cLabel;
                        %need to add in multiplicities
                        
                        SwapData(5,SwapData(4,:)==dLabel) = SwapBackward;
                        SwapData(4,SwapData(4,:)==dLabel) = aLabel;
                        %need to add in multiplicities
                        
                        MultLabel = dLabel;
                        
                        AdjustPath = [AdjustPath;[{'pDownInv';'RCW';'RCCWInv';'Swap'},{[cLabel;-cDir];...
                            [bLabel, cLabel, aLabel, MultLabel,MultLabel; bDir,cDir,-aDir,0,0];...
                            [aLabel, bLabel, cLabel, MultLabel,MultLabel; aDir,bDir,-cDir,0,0];SwapData}]];
                        %disp('Adjusted')
                        
                    end
                end
                
                %{
                for ll = (length(Path)-1):-1:1
                    %Note that Path(end) = 0
                     
                    if IsTopPath(ll)
                        Loc = [1,2*Path(ll)+1, 2*Path(ll)+2,2];
                        Loc2 = Loc([2,3,4,1]);
                        
                        dLabel = StructureTemp(Loc(1));
                        aLabel = StructureTemp(Loc(2));
                        bLabel = StructureTemp(Loc(3));
                        cLabel = StructureTemp(Loc(4));
                        
                        aDir = ChargeDirectionsTemp(Loc(2));
                        bDir = ChargeDirectionsTemp(Loc(3));
                        cDir = ChargeDirectionsTemp(Loc(4));
                        dDir = -aDir;
                        
                        StructureTemp(Loc) = StructureTemp(Loc2);
                        ChargeDirectionsTemp(Loc(1)) = -ChargeDirectionsTemp(Loc(2));
                        ChargeDirectionsTemp(Loc) = ChargeDirectionsTemp(Loc2);
                        
                        SwapData = SwapDataOriginal;
                        SwapData(2,SwapData(1,:)==dLabel) = 1;
                        SwapData(1,SwapData(1,:)==dLabel) = aLabel;
                        %need to add in multiplicities
                        
                        SwapData(5,SwapData(4,:)==dLabel) = 1;
                        SwapData(4,SwapData(4,:)==dLabel) = cLabel;
                        %need to add in multiplicities
                        
                        MultLabel = dLabel;
                        
                        AdjustPath = [AdjustPath;[{'Swap';'RCCW';'RCWInv';'pDown'},{SwapData;[dLabel, bLabel, cLabel, MultLabel,MultLabel; -dDir,bDir,-cDir,0,0];...
                            [bLabel, cLabel, dLabel, MultLabel,MultLabel; bDir,cDir,dDir,0,0];[cLabel;-cDir]}]];
                        %disp('AdjustedOther2')
                        
                    else
                        Loc = [1,2*Path(ll)+1, 2*Path(ll)+2,2];
                        Loc2 = Loc([4,1,2,3]);
                        
                        dLabel = StructureTemp(Loc(4));
                        aLabel = StructureTemp(Loc(1));
                        bLabel = StructureTemp(Loc(2));
                        cLabel = StructureTemp(Loc(3));
                        
                        aDir = ChargeDirectionsTemp(Loc(1));
                        bDir = ChargeDirectionsTemp(Loc(2));
                        cDir = ChargeDirectionsTemp(Loc(3));
                        dDir = ChargeDirectionsTemp(Loc(4));%-cDir;
                        
                        StructureTemp(Loc) = StructureTemp(Loc2);
                        ChargeDirectionsTemp(Loc(4)) = -ChargeDirectionsTemp(Loc(3));
                        ChargeDirectionsTemp(Loc) = ChargeDirectionsTemp(Loc2);
                        
                        SwapData = SwapDataOriginal;
                        SwapData(2,SwapData(1,:)==dLabel) = 1;
                        SwapData(1,SwapData(1,:)==dLabel) = cLabel;
                        %need to add in multiplicities
                        
                        SwapData(5,SwapData(4,:)==dLabel) = 1;
                        SwapData(4,SwapData(4,:)==dLabel) = aLabel;
                        %need to add in multiplicities
                        
                        MultLabel = dLabel;
                        
                        AdjustPath = [AdjustPath;[{'pDownInv';'RCW';'RCCWInv';'Swap'},{[cLabel;-cDir];...
                            [bLabel, cLabel, dLabel, MultLabel,MultLabel; bDir,cDir,dDir,0,0];...
                            [dLabel, bLabel, cLabel, MultLabel,MultLabel; -dDir,bDir,-cDir,0,0];SwapData}]];
%                        AdjustPath = [AdjustPath;[{'pDownInv';'pDown';'RCW';'RCCWInv';'Swap'},{[cLabel;-cDir];[aLabel;aDir];...
%                            [bLabel, cLabel, dLabel, MultLabel,MultLabel; bDir,cDir,dDir,0,0];...
%                            [dLabel, bLabel, cLabel, MultLabel,MultLabel; -dDir,bDir,-cDir,0,0];SwapData}]];
                        %disp('Adjusted2')
                        
                    end
                end
                %}
                RotationModifications{kk} = [RotationModifications{kk};AdjustPath];
                
                
                %now want to adjust so the internal bonds are counting up
                %from left to right
                
                if numel(StructureTemp) > 2
                    WorkFrom = StructureTemp(1,StructureTemp(2,:) == -LegsOutOrder(end));
                    Ordering = zeros([1,0]);
                    
                    %fix this
                    
                    while ~isempty(WorkFrom)
                        Current = WorkFrom(end);
                        Root = find(any(StructureTemp == WorkFrom(end),1))-1;
                        if numel(Root)~=1
                            error('Affirmation Error: Root should be exactly one entry')
                        end
                        WorkFrom(end) = [];
                        
                        %{
                        if any(Root == StructureTemp(2,:))
                            Ordering = [Ordering,Current];
                        else
                            Ordering = [Current,Ordering];
                        end
                        %}
                        Ordering = [Current,Ordering];
                        
                        End = StructureTemp(:,Current+1)';
                        WorkFrom = [WorkFrom,End(End>0)];
                    end
                    
                    if ~isequal(sort(Ordering,'ascend'), 1:numel(ChargeSideInt))
                        error('Affirmation Error: ordering should be equal to the number of internal entries')
                    end
                    
                    %if this is different from ideal then:
                    if ~isequal(Ordering, 1:numel(Ordering))
                        
                        SwapData = SwapDataOriginal;
                        [~,InverseOrdering] = sort(Ordering,'ascend');
                        SwapData(1,(numel(ChargeSide)+1):end) = Ordering;
                        SwapData(3,1:numel(ChargeSideInt)) = Ordering;
                        SwapData(4,(numel(ChargeSide)+1):end) = InverseOrdering;
                        SwapData(6,1:numel(ChargeSideInt)) = InverseOrdering;
                        
                        StructureNew = StructureTemp;
                        for pp = 1:numel(Ordering)
                            StructureNew(StructureTemp == Ordering(pp)) = pp;
                        end
                        
                        StructureTemp = StructureNew(:,[1,Ordering+1]);
                        ChargeDirectionsTemp = ChargeDirectionsTemp(:,[1,Ordering+1]);
                        
                        RotationModifications{kk} = [RotationModifications{kk};{'Swap',SwapData}];
                        
                    end
                    
                end
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %fifth part is to use the F-moves to get the points to be most
                %right
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                for ff = (size(LegsOutOrder,2)-1):-1:3
                    %now we want to make sure the kk-2 vertex is [kk-3,-kk]
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %now work out F-moves
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    LocUsing = find(any(StructureTemp(2,:)==(-LegsOutOrder(ff)),1))-1; %this will throw an error if it is not on the bottom (i.e. right most)
                    if isempty(LocUsing)
                        error('Affirmation Error: we should be able to find at least one vertex which should be made right-%i th-most',ff)
                    end
                    if numel(LocUsing)>1
                        error('Affirmation Error: there should only be one vertex which should be made right-%i th-most',ff)
                    end
                    
                    
                    Path = LocUsing;
                    %don't need IsTopPath as it should always be on the right
                   
                    while LocUsing ~= ff-1
                        LocUsing = find(any(StructureTemp == LocUsing,1))-1;
                        Path = [Path, LocUsing];
                        if isempty(LocUsing)
                            error('Affirmation Error: We couldn''t find the previous vertex from this Loc Using')
                        elseif numel(LocUsing)>1
                            error('Affirmation Error: We found too many previous verticies, we should only have found one')
                        end
                    end
                    
                    
                    for ll = (length(Path)-1):-1:1
                        
                        %Note that Path(end) = kk-3
                        aLabel = StructureTemp(1,Path(end)+1);
                        bLabel = StructureTemp(1,Path(ll)+1);
                        cLabel = StructureTemp(2,Path(ll)+1);
                        dLabel = Path(end);
                        eLabel = Path(ll);
                        
                        %verify
                        if StructureTemp(2,Path(end)+1) ~= eLabel
                            error('Affirmation Error: the path seems to be wrong in this calculation.')
                        end
                        
                        muMult = Path(ll);
                        nuMult = Path(end);
                        
                        aDirection = ChargeDirectionsTemp(1,Path(ll+1)+1);
                        bDirection = ChargeDirectionsTemp(1,Path(ll)+1);
                        cDirection = ChargeDirectionsTemp(2,Path(ll)+1);
                        
                        %d and e are going to be internal legs so they are
                        %already set to 1
                        
                        %adjust the StructureTemp, and ChargeDirectionsTemp
                        StructureTemp(2,Path(end)+1) = StructureTemp(2,Path(ll)+1);
                        StructureTemp(2,Path(ll)+1) = StructureTemp(1,Path(ll)+1);
                        StructureTemp(1,Path(ll)+1) = StructureTemp(1,Path(end)+1);
                        StructureTemp(1,Path(end)+1) = Path(ll);
                        
                        ChargeDirectionsTemp(2,Path(end)+1) = ChargeDirectionsTemp(2,Path(ll)+1);
                        ChargeDirectionsTemp(2,Path(ll)+1) = ChargeDirectionsTemp(1,Path(ll)+1);
                        ChargeDirectionsTemp(1,Path(ll)+1) = ChargeDirectionsTemp(1,Path(end)+1);
                        ChargeDirectionsTemp(2,Path(end)+1) = +1; %this is f which is an internal leg in all cases (so its going up by our fixed convention).
                        
                        %add the F-move:
                        RotationModifications{kk} = [RotationModifications{kk}; {'FInv',...
                            [aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,muMult,nuMult,muMult,nuMult;aDirection,bDirection,cDirection,+1,+1,+1,0,0,0,0]}];
                        
                    end
                    
                    NewStructurePermute = reshape(1:numel(StructureTemp), size(StructureTemp));
                    
                    Have = StructureTemp(1,ff-1+1);
                    if Have ~= -LegsOutOrder(ff)
                        if Have<0
                            error('Affirmation Error: We were suppose to find an internal leg here or the wanted output leg, this is for external leg %i',ff)
                        end
                        Wanted = ff-2;
                        
                        if Have ~= Wanted
                            HaveLocation = find(StructureTemp == Have);
                            WantedLocation = find(StructureTemp == Wanted);
                            
                            if isempty(HaveLocation)
                                error('Affirmation Error: can''t find the Have location for the %i th external leg',ff)
                            end
                            if numel(HaveLocation)>1
                                error('Affirmation Error: can find more then one entry for the Have location for the %i th external leg',ff)
                            end
                            
                            
                            if isempty(WantedLocation)
                                error('Affirmation Error: can''t find the Have location for the %i th external leg',ff)
                            end
                            if numel(WantedLocation)>1
                                error('Affirmation Error: can only find one entry for the Have location for the %i th external leg',ff)
                            end
                            
                            
                            
                            if Have<0 || Wanted <0
                                error('Affirmation Error: Have and Wanted should be internal');
                            end
                            
                            StructureTemp(HaveLocation) = Wanted;
                            StructureTemp(WantedLocation) = Have;
                            StructureTemp(:,[Have,Wanted]+1) = StructureTemp(:,[Wanted,Have]+1);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %now do the swap
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            NewStructurePermute = sort(StructureTemp(:)','ascend');
                            PosNumbers = sum(NewStructurePermute>0);
                            NegNumbers = sum(NewStructurePermute<0);
                            NewStructurePermute(1:NegNumbers) = NewStructurePermute(NegNumbers:-1:1);
        
                            NewStructurePermute(NegNumbers+[Have,Wanted]) = [Wanted,Have];
                            NewStructurePermuteMult = [NewStructurePermute((NegNumbers+1):end),zeros([1,NegNumbers])];
                            
                            NewStructureDual = false(size(NewStructurePermute)); %everything is +1
                            
                            
                            
                            DetailsSwap = [NewStructurePermute; NewStructureDual; NewStructurePermuteMult;...
                                   NewStructurePermute; NewStructureDual; NewStructurePermuteMult];
            
                            RotationModifications{kk} = [RotationModifications{kk}; {'Swap',DetailsSwap}];
                        end     
                    else
                        error('Affirmation Error: We found the %i th external leg was on the left when it should have been right most',ff)
                        
                        %Have = StructureTemp(2,ff-1+1);
                        %Wanted = ff-2;
                    end
                    
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %sixth part is to relabel external charges and store final
                %data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                EntriesList = [StructureTemp(1,2),StructureTemp(2,2)];
                EntriesList = [EntriesList, StructureTemp(2,3:end),StructureTemp(2,1)];
                    
                
                if 2*numel(EntriesList)-2~=numel(StructureTemp)
                    error('Affirmation Error: We seem to have the wrong number of Entries in EntriesList');
                end
                
                if any(EntriesList>0)
                    error('Affirmation Error: The EntriesList seem to have internal legs when they should never have this');
                end
                
                if ~isequal(sort(-EntriesList,'ascend'), 1:numel(EntriesList))
                    error('Affirmation Error: The EntriesList does not seem to have all the external legs');
                end
                
                
                if ~isequal(EntriesList, -(1:numel(EntriesList)))
                    NewStructurePermute = sort(StructureTemp(:)','ascend');
                    PosNumbers = sum(NewStructurePermute>0);
                    NegNumbers = sum(NewStructurePermute<0);
                    NewStructurePermuteInv = NewStructurePermute;
                    
                    InverseLabelingIndex = -EntriesList;
                    [~,LabelingIndex] = sort(-EntriesList,'ascend');
                    
                    NewStructurePermuteInv(1:NegNumbers) = -LabelingIndex;
                    NewStructurePermute(1:NegNumbers) = -InverseLabelingIndex;
                    
                    NewStructureDual = false(size(NewStructurePermute)); %everything is +1
                    NewStructurePermuteMult = [NewStructurePermute((NegNumbers+1):end),zeros([1,NegNumbers])];
                    
                    
                    DetailsSwap = [NewStructurePermute; NewStructureDual; NewStructurePermuteMult;...
                           NewStructurePermuteInv; NewStructureDual; NewStructurePermuteMult];
                    
                    RotationModifications{kk} = [RotationModifications{kk}; {'Swap',DetailsSwap}];
                end
                
                end
            end
            
            %Last part is to invert if required
            
            if Invert
                
                for kk = 1:numel(RotationModifications)
                    NewRotationModifications = cell(size(RotationModifications{kk}));
                    OldRotationModifications = RotationModifications{kk};
                    
                    for ll = 1:size(RotationModifications{kk},1)
                        
                        switch RotationModifications{kk}{ll,1}
                            
                            case 'F'
                                NewRotationModifications{ll,1} = 'FInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                                
                            case 'FInv'
                                NewRotationModifications{ll,1} = 'F';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'FSq'
                                NewRotationModifications{ll,1} = 'FSqInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'FSqInv'
                                NewRotationModifications{ll,1} = 'FSq';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                                
                                
                                
                                
                            case 'R'
                                NewRotationModifications{ll,1} = 'RInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                %NewRotationModifications{ll,2}(:,[1,2]) = NewRotationModifications{ll,2}(:,[2,1]);
                                
                            case 'RInv'
                                NewRotationModifications{ll,1} = 'R';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                %NewRotationModifications{ll,2}(:,[1,2]) = NewRotationModifications{ll,2}(:,[2,1]);
                                
                                
                            case 'B' %these shouldn't show up but are here for completeness (if I want to reuse this code later)
                                NewRotationModifications{ll,1} = 'BInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                %NewRotationModifications{ll,2}(:,[2,3]) = NewRotationModifications{ll,2}(:,[3,2]);
                                
                                
                            case 'BInv' %these shouldn't show up but are here for completeness (if I want to reuse this code later)
                                NewRotationModifications{ll,1} = 'B';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                %NewRotationModifications{ll,2}(:,[2,3]) = NewRotationModifications{ll,2}(:,[3,2]);
                                
                                
                                
                                
                            case 'RCW'
                                NewRotationModifications{ll,1} = 'RCWInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'RCWInv'
                                NewRotationModifications{ll,1} = 'RCW';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'RCCW'
                                NewRotationModifications{ll,1} = 'RCCWInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'RCCWInv'
                                NewRotationModifications{ll,1} = 'RCCW';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'pDown'
                                NewRotationModifications{ll,1} = 'pDownInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'pDownInv'
                                NewRotationModifications{ll,1} = 'pDown';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                            
                            case 'pUp'
                                NewRotationModifications{ll,1} = 'pUpInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'pUpInv'
                                NewRotationModifications{ll,1} = 'pUp';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'Round'
                                NewRotationModifications{ll,1} = 'RoundInv';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'RoundInv'
                                NewRotationModifications{ll,1} = 'Round';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                            case 'Swap'
                                %the swap funciton 
                                NewRotationModifications{ll,1} = 'Swap';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2}([4:6,1:3],:);
                                
                                
                            case 'FlipChargeDirections'
                                NewRotationModifications{ll,1} = 'FlipChargeDirections';
                                NewRotationModifications{ll,2} = OldRotationModifications{ll,2};
                                
                                
                                
                            case 'Dim' %these shouldn't show up but are here for completeness (if I want to reuse this code later)
                                error('Affirmation Error: I can''t invert the Dim case')
                                
                            otherwise
                                error('Affirmation Error: We couldn''t find one of the entries')
                            
                            
                        end
                        
                        
                    end
                    
                    RotationModifications{kk} = NewRotationModifications(end:-1:1,:);
                    
                end
                
                
            end
            
            if FlagSingleOutput
                RotationModifications = RotationModifications{1};
            end
            
        end
        
        function Structure = GenerateStructure(Legs)
            
            TotalLegsBottom = Legs(1);
            TotalLegsTop = Legs(2);
            
            if TotalLegsBottom*TotalLegsTop >0
                Structure = zeros([2,TotalLegsBottom+TotalLegsTop-1]);
                
                if TotalLegsBottom == 1
                    Structure(1,1) = -1;
                else
                    Structure(1,1) = TotalLegsBottom-1;
                    Structure(2,TotalLegsBottom) = -TotalLegsBottom;
                end
                
                for jj = 2:(TotalLegsBottom-1)
                    Structure(1,TotalLegsBottom-jj+2) = TotalLegsBottom-jj;
                    Structure(2,TotalLegsBottom-jj+1) = -(TotalLegsBottom-jj+1);
                end
                
                if TotalLegsBottom>=2
                    Structure(1,2) = -1;
                end
                
                if TotalLegsTop == 1
                    Structure(2,1) = -TotalLegsBottom-1;
                else
                    Structure(2,1) = TotalLegsBottom+TotalLegsTop-2;
                    Structure(2,TotalLegsBottom+TotalLegsTop-1) = -(TotalLegsBottom+TotalLegsTop);
                end
                
                for jj = 2:(TotalLegsTop-1)
                    Structure(1,TotalLegsTop-jj+2+TotalLegsBottom-1) = TotalLegsTop-jj-1+TotalLegsBottom;
                    Structure(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = -(TotalLegsTop-jj+1+TotalLegsBottom);
                end
                
                if TotalLegsTop>=2
                    Structure(1,TotalLegsBottom+1) = -(TotalLegsBottom+1);
                end
                
            else
                
                if TotalLegsBottom == 0
                    
                    if TotalLegsTop == 2
                        Structure(1,1) = -1;
                        Structure(2,1) = -2;
                    else
                        Structure(1,1) = TotalLegsTop-2;
                        Structure(2,1) = -TotalLegsTop;
                        Structure(1,2) = -1;
                        Structure(2,2) = -2;
                    end
                     
                    for jj = 2:(TotalLegsTop-2)
                        Structure(1,TotalLegsTop-jj+1) = TotalLegsTop-jj-1;
                        Structure(2,TotalLegsTop-jj+1) = -(TotalLegsTop-jj+1);
                    end
                    
                else %if TotalLegsTop == 0
                    
                    if TotalLegsBottom == 2
                        Structure(1,1) = -1;
                        Structure(2,1) = -2;
                    else
                        Structure(1,1) = TotalLegsBottom-2;
                        Structure(2,1) = -TotalLegsBottom;
                        Structure(1,2) = -1;
                        Structure(2,2) = -2;
                    end
                     
                    for jj = 2:(TotalLegsBottom-2)
                        Structure(1,TotalLegsBottom-jj+1) = TotalLegsBottom-jj-1;
                        Structure(2,TotalLegsBottom-jj+1) = -(TotalLegsBottom-jj+1);
                    end
                    
                end
            
            end
            
            if TotalLegsTop+TotalLegsBottom<2
                Structure = zeros([2,0]);
            end
            
        end
        
        function StoredLocations = GenerateStoredLocations(Structure)
            if max(Structure(:))>0
                
                for jj = max(Structure(:)):-1:1
                    Updating = jj;
                    Using = [];
                    while ~isempty(Updating)
                        kk = Updating(end);
                        Updating(end) = [];
                        Use1 = Structure(1,kk+1);
                        Use2 = Structure(2,kk+1);
                        if Use1>0
                            if Use1>jj
                                Using = [Using,StoredLocations{Use1}];
                            else
                                Updating = [Updating,Use1];
                            end
                        else %Use1<0
                            Using = [Using, -Use1];
                        end
                        
                        if Use2>0
                            if Use2>jj
                                Using = [Using,StoredLocations{Use2}];
                            else
                                Updating = [Updating,Use2];
                            end
                        else %Use2<0
                            Using = [Using, -Use2];
                        end
                    end
                    StoredLocations{jj} = Using;
                end
            else
                StoredLocations = cell([0,1]);
            end
        end
        
    end
    
    methods(Access = 'public')
        
        function [Eigenvalue, Eigenvector] = MaxEigenValue(self, Legs, Shape)
            
            if nargin<2
                Legs = [];
            end
            
            FlagAbs = false;
            
            %if ~self.SymHandle.IsSelfDual
                InverseIrrep = self.SymHandle.getInverseIrrep;
                if isequal(self.ChargeSide, [-ones([1,Legs(1)]),ones([1,Legs(2)])])
                    FlagCorrectSides = true;
                elseif isequal(self.ChargeSide, [ones([1,Legs(1)]),-ones([1,Legs(2)])])
                    FlagCorrectSides = false;
                    self.ChargeSide = -self.ChargeSide;
                    self.ChargeLabelsInternal = InverseIrrep(self.ChargeLabelsInternal);
                    self.ChargeDirections(self.Structure>0) = -self.ChargeDirections(self.Structure>0);
                else
                    error('Error: The Tensor MaxEigenValue doesn''t work unless all the legs of a side are grouped up at any one time')
                end
            %end
            if nargin<3
                [Matrix,Data] = SymTensor.SymTen2Mat(self, Legs);
            else
                [Matrix,Data] = SymTensor.SymTen2Mat(self, Legs,Shape);
            end
            
            if ~isempty(Data{8})
                error('Error: Tensor MaxEigenValue doesn''t work with currently braided tensors, the unresolved braiding has to be removed before we can proceed');
            end
            
            Flavours = cell([1,size(Matrix,1)]);
            CellUsed = cell([1,size(Matrix,1)]);
            MaxEigenvector = cell([1,size(Matrix,1)]);
            
            Numbers = ([1,size(Matrix,1)]);
            MaxEigenvalues = zeros(size(Matrix));
            for kk = 1:numel(Matrix)
                 if size(Matrix{kk},1)~=size(Matrix{kk},2)
                     error('Error: this is not a square Tensor')
                 end
                 
                 if FlagAbs
                     Eigenvalues = abs(eig(Matrix{kk}));
                 else
                     Eigenvalues = real(eig(Matrix{kk}));
                 end
                 MaxEigenvalues(kk) = max(Eigenvalues);
                 MaxEigenvector{kk} = mat2cell(Matrix{kk}(:,MaxEigenvalues(kk)==Eigenvalues),size(Matrix{kk},1), ones([1,sum(MaxEigenvalues(kk)==Eigenvalues)]));
                 Numbers(kk) = sum(MaxEigenvalues(kk)==Eigenvalues);
                 CellUsed{kk} = repmat(kk,[1,Numbers(kk)]);
                 Flavours{kk} = repmat(Data{1}(kk),[1,Numbers(kk)]);
            end
            clear Matrix
            Eigenvalue = max(MaxEigenvalues);
            %CellUsed = 1:numel(Data{1});
            %Flavours = Data{1};
            
            CellUsed = CellUsed(MaxEigenvalues == Eigenvalue);
            MaxEigenvector = MaxEigenvector(MaxEigenvalues == Eigenvalue);
            Flavours = Flavours(MaxEigenvalues == Eigenvalue);
            
            %check that we are getting only a single vector
            ChargesExt = Data{7}(Data{11}<0);
            [~,Index] = sort(-Data{11}(Data{11}<0),'ascend');
            ChargesExt = ChargesExt(Index);
            
            
            %recreate the following
            Data{3} = {1};
            Data{4}(2) = 1;
            Data{7} = -ones([2,Data{4}(1)]);
            Data{7}(2,1) = +1;
            if Data{4}(1)>1
                Data{7}(1,2) = ChargesExt(1);
                Data{7}(2,2:end) = ChargesExt(2:Data{4}(1));
            else
                Data{7}(1,1) = ChargesExt(1);
            end
            Data{11} = SymTensor.GenerateStructure(Data{4});
            
            
            %Now duplicate
            MaxEigenvector = horzcat(MaxEigenvector{:});
            Flavours = horzcat(Flavours{:});
            CellUsed = horzcat(CellUsed{:});
            
            
            Eigenvector = cell(size(MaxEigenvector));
            %The following we select out something in particula
            DataSet = repmat({Data},size(MaxEigenvector));
            for kk = 1:numel(DataSet)
                DataSet{kk}{1} = Flavours(kk);
                DataSet{kk}{2} = DataSet{kk}{2}(CellUsed(kk));
                DataSet{kk}{5} = DataSet{kk}{5}(CellUsed(kk));
                DataSet{kk}{6} = DataSet{kk}(1);
                
                Eigenvector{kk} = SymTensor.SymMat2Ten(MaxEigenvector(kk),DataSet{kk});
            end
            
            
            if numel(MaxEigenvector) == 1
                Eigenvector = Eigenvector{1};
            end
            
            %Eigenvector = 0;
        end
        
        function self = ForceTrivial(self,Legs, RemoveExternalLegs)
            %this function 
            
            if nargin<3
                RemoveExternalLegs = true;
            end
            
            LegsInternal = Legs(Legs>0);
            LegsExternal = -Legs(Legs<0);
            TrivialNumber = self.SymHandle.TrivialIrrep;
            RemovedLegsExternal = LegsExternal;
            
            %first fix all internal legs we want
            Keep = all(self.ChargeLabelsInternal(:,LegsInternal) == TrivialNumber,2);
            self.TensorEntries = self.TensorEntries(Keep);
            self.TensorEntriesSizes = self.TensorEntriesSizes(Keep,:);
            self.ChargeLabelsExternal = self.ChargeLabelsExternal(Keep,:);
            self.ChargeLabelsInternal = self.ChargeLabelsInternal(Keep,:);
            self.MultiplicitiesInternal = self.MultiplicitiesInternal(Keep,:);
            
            %now adjust the external legs
            Keep = all(self.ChargeLabelsExternal(:,LegsExternal) == TrivialNumber,2);
            self.TensorEntries = self.TensorEntries(Keep);
            self.TensorEntriesSizes = self.TensorEntriesSizes(Keep,:);
            self.ChargeLabelsExternal = self.ChargeLabelsExternal(Keep,:);
            self.ChargeLabelsInternal = self.ChargeLabelsInternal(Keep,:);
            self.MultiplicitiesInternal = self.MultiplicitiesInternal(Keep,:);
            
            if RemoveExternalLegs && ~isempty(RemovedLegsExternal)
                
                
                RemovedLegsInternal = zeros([1,0]);
                ExternalReplacement = zeros([1,0]);
                NumberStored = zeros([1,0]);
                ChargeValue = zeros([1,0]);
                for kk = 1:numel(self.StoredLocations)
                    Counts = any(repmat(self.StoredLocations{kk}(:),[1,numel(RemovedLegsExternal)]) == repmat(RemovedLegsExternal(:)',[numel(self.StoredLocations{kk}),1]),2);
                    if sum(Counts,1)==numel(self.StoredLocations{kk})
                        RemovedLegsInternal = [RemovedLegsInternal,kk-1];
                        ExternalReplacement = [ExternalReplacement,0];
                        NumberStored = [NumberStored,numel(self.StoredLocations{kk})];
                        ChargeValue = [ChargeValue,0];
                    elseif sum(Counts,1)==numel(self.StoredLocations{kk})-1
                        RemovedLegsInternal = [RemovedLegsInternal,-(kk-1)];
                        ExternalReplacement = [ExternalReplacement,self.StoredLocations{kk}(~Counts)];
                        NumberStored = [NumberStored,numel(self.StoredLocations{kk})];
                        NewCharge = self.ChargeDirections(self.Structure == -self.StoredLocations{kk}(~Counts));
                        ChargeValue = [ChargeValue,NewCharge];
                    end
                end
                
                
                
                
                Temp = self.CyclicLegs;
                Temp(Temp == RemovedLegsExternal) = [];
                [~,self.CyclicLegs] = sort(Temp,'ascend');
                [~,self.CyclicLegs] = sort(self.CyclicLegs,'ascend');
                
                
                Temp = self.CyclicLegs;
                Temp(Temp == RemovedLegsExternal) = [];
                [Temp,self.CyclicLegsBraided] = sort(Temp,'ascend');
                [~,self.CyclicLegsBraided] = sort(self.CyclicLegsBraided,'ascend');
                
                [~,NewName] = sort([Temp,RemovedLegsExternal],'ascend');
                [~,NewName] = sort(NewName,'ascend');
                
                
                if size(self.ChargeLabelsExternal,2)-numel(RemovedLegsExternal)<=1
                    %if we have 1 or less external leg out we delete
                    %everything (we will have no InternalLegs left)
                    self.Structure = [];
                    self.ChargeDirections = [];
                else
                    %otherwise we need to consider the core, unlike the
                    %internal legs we can't remove the core if we have more
                    %then 1 leg
                    
                    %if we removed an internal leg off with a replacement
                    %then we aren't really removing the core, otherwise we are 
                    
                    TopRemoved = any(self.Structure(1,1) == [RemovedLegsInternal,-RemovedLegsExternal],1);
                    TopReplaced = ~TopRemoved;
                    Replaced = TopReplaced||any(self.Structure(2,1) == [RemovedLegsInternal,-RemovedLegsExternal],1);
                    
                    
                    if (sum(self.Structure(1,1) == [RemovedLegsInternal,-RemovedLegsExternal])+sum(self.Structure(1,1) == [RemovedLegsInternal,-RemovedLegsExternal]))>1
                        error('Affirmation Error: the core is having issues')
                    end
                    
                    if Replaced && ~(any(self.Structure(2-TopReplaced,1) == RemovedLegsInternal) &&...
                            ExternalReplacement(self.Structure(2-TopReplaced,1) == RemovedLegsInternal)~=0)
                        %then everything is not safe
                        
                        %the other one should be internal, move it in place
                        %of the core and don't provide a replacement
                        
                        InternalLegToCore = self.Structure(2-TopReplaced,1);
                        if InternalLegToCore <0
                            error('Affirmation Error: We shouldn''t have gotten an internal leg');
                        end
                        
                        self.Structure(:,1) = self.Structure(:,InternalLegToCore+1);
                        self.ChargeDirections(:,1) = self.ChargeDirections(:,InternalLegToCore+1);
                        
                        
                        RemovedLegsInternal = [RemovedLegsInternal,InternalLegToCore];
                        ExternalReplacement = [ExternalReplacement,0];
                        NumberStored = [NumberStored,0];
                        ChargeValue = [ChargeValue,0];
                        
                        self.Structure(:,RemovedLegsInternal+1) = [];
                        self.ChargeDirections(:,RemovedLegsInternal+1) = [];
                        
                    end
                    
                    %{
                    if (RemoveLegsInternal~=0)
                        %now we know the core hase been removed and needs
                        %to be replaced
                        FlagRemoveCore = true;
                        TempCore = self.Structure(:,1)(all(self.Structure(:,1) ~= [RemoveLegsInternal,-RemoveLegsExternal]));
                        CountsCore = any(repmat(self.StoredLocations{kk}(:),[1,numel(RemovedLegsExternal)]) == repmat(RemovedLegsExternal(:)',[numel(self.StoredLocations{kk}),1]),2);
                        RemovedLegsInternal = [RemovedLegsInternal,TempCore];
                        ExternalReplacement = [ExternalReplacement,self.StoredLocations{TempCore}(~CountsCore)];
                        %NumberStored = [NumberStored,numel(self.StoredLocations{kk})];
                        %NewCharge = self.ChargeDirections(self.Structure == -self.StoredLocations{kk}(~Counts));
                        ChargeValue = [ChargeValue,NewCharge];
                        
                    else
                        FlagRemoveCore = false;
                        
                    end
                    %}
                    
                    self.TensorEntriesSizes(:,RemovedLegsExternal) = [];
                    self.ChargeLegDimensions(:,RemovedLegsExternal) = [];
                    self.ChargeSide(RemovedLegsExternal) = [];
                    self.ChargeSideInt(RemovedLegsInternal) = [];
                    
                    %replacements for Structure and ChargeDirections
                    Re_add = unique(ExternalReplacement);
                    
                    for kk = ExternalReplacement
                        if kk ~=0
                            MaxNumberStored = max(NumberStored(RemovedLegsInternal == kk));
                            ListNumber = find((NumberStored == MaxNumberStored) & (RemovedLegsInternal == kk));
                            
                            self.ChargeDirections(self.Structure == RemovedLegsInternal(ListNumber)) = ChargeValue(RemovedLegsInternal);
                            self.Structure(self.Structure == RemovedLegsInternal(ListNumber)) = -ExternalReplacement(RemovedLegsInternal);
                        end
                    end
                    
                    self.StoredLocations(RemovedLegsInternal) = [];
                    for kk = 1:numel(self.StoredLocations)
                        self.StoredLocations{kk} = NewName(self.StoredLocations{kk});
                    end
                    
                    
                    %need newname
                    %need to fix the external replacement section
                    %need to modify structure to pick out whenever one of
                    %its inputs have been removed
                    
                    
                    
                end
                
                self.ChargeLabelsExternal(:,RemovedLegsExternal) = [];
                self.ChargeLabelsInternal(:,RemovedLegsInternal) = [];
                self.MultiplicitiesInternal(:,RemovedLegsInternal) = [];
                
                %self.Braiding
                %self.BraidingDirection
                
            end
            
            %need to add rotations for internal vetricies that get mixed
            %up(this will also account for the change in normalisation)
            
        end
        
        function [U,Sig,V] = SVD(self, Legs, Shape)
            if nargin<2
                Legs = [];
            end
            
            %if ~self.SymHandle.IsSelfDual
                InverseIrrep = self.SymHandle.getInverseIrrep;
                if isequal(self.ChargeSide, [-ones([1,Legs(1)]),ones([1,Legs(2)])])
                    FlagCorrectSides = true;
                elseif isequal(self.ChargeSide, [ones([1,Legs(1)]),-ones([1,Legs(2)])])
                    FlagCorrectSides = false;
                    self.ChargeSide = -self.ChargeSide;
                    self.ChargeLabelsInternal = InverseIrrep(self.ChargeLabelsInternal);
                    self.ChargeDirections(self.Structure>0) = -self.ChargeDirections(self.Structure>0);
                else
                    error('Error: The Tensor SVD doesn''t work unless all the legs of a side are grouped up at any one time')
                end
            %end
            if nargin<3
                [Matrix,Data] = SymTensor.SymTen2Mat(self, Legs);
            else
                [Matrix,Data] = SymTensor.SymTen2Mat(self, Legs,Shape);
            end
            
            if ~isempty(Data{8})
                error('Error: Tensor SVD doesn''t work with currently braided tensors, the unresolved braiding has to be removed before we can proceed');
            end
            
            UMat = cell(size(Matrix));
            VMat = cell(size(Matrix));
            SigMat = cell(size(Matrix));
            
            DimensionMat{1} = cell(size(Data{1}));
            DimensionMat{2} = DimensionMat{1};
            
            for nn = 1:numel(Matrix)
                [UMat{nn}, SigMat{nn}, VMat{nn}] = svd(Matrix{nn});
                DimensionMat1{nn} = size(Matrix{nn},1);
                DimensionMat2{nn} = size(Matrix{nn},2);
            end
            
            
            ChargeDirectionsTotal = Data{7};
            StructureTotal = Data{11};
            
            SigChargeDirections = [-1;1];
            SigStructure = [-1;-2];
            
            UStructure = zeros([2,Legs(1)]);
            UChargeDirections = -ones([2,Legs(1)]);
            UChargeDirections(2,1) = 1;
            UStructure(2,1) = -(Legs(1)+1);
            for jj = 1:(Legs(1)-2)
                UStructure(1,Legs(1)-jj+1) = jj;
                UStructure(2,Legs(1)-jj+1) = -Legs(1)+jj-1;
                UChargeDirections(2,Legs(1)-jj+1) = ChargeDirectionsTotal(StructureTotal == (-Legs(1)+jj-1));
            end
            if Legs(1) > 1
                UStructure(1,1) = Legs(1)-1;
                UStructure(1,2) = -1;
                UStructure(2,2) = -2;
                
            elseif Legs(1) == 1
                UStructure(1,1) = -1;
            else%if Legs(1) == 0
                error('Not doing this')
            end
            
            VStructure = zeros([2,Legs(2)]);
            VChargeDirections = -ones([2,Legs(2)]);
            VChargeDirections(2,1) = 1;
            VStructure(2,1) = -(Legs(2)+1);
            for jj = 1:(Legs(2)-2)
                VStructure(1,Legs(2)-jj+1) = jj;
                VStructure(2,Legs(2)-jj+1) = -Legs(2)+jj-1;
                VChargeDirections(2,Legs(2)-jj+1) = -ChargeDirectionsTotal(StructureTotal == (-Legs(2)+jj-1   -Legs(1)));
            end
            if Legs(2) > 1
                VStructure(1,1) = Legs(2)-1;
                VStructure(1,2) = -1;
                VStructure(2,2) = -2;
                
            elseif Legs(2) == 1
                VStructure(1,1) = -1;
            else%if Legs(2) == 0
                error('Not doing this')
            end
            
            UUniqueEntries = Data{1};
            VUniqueEntries = Data{1};
            
            
            DataInternal = mat2cell(Data{1},ones(size(Data{1})),1);
            
            DataU = {UUniqueEntries, Data{2},DimensionMat1, [Legs(1),1], Data{5}, DataInternal,...
                UChargeDirections, [],[],Data{10},UStructure,0};%, [], [], , Data{11}, , };
            DataV = {VUniqueEntries, Data{3},DimensionMat2, [Legs(2),1], Data{6}, DataInternal,...
                VChargeDirections, [],[],Data{10},VStructure,0};%, [], [], , Data{11}, , };
            DataSig = {Data{1}, DimensionMat1, DimensionMat2, [1,1], DataInternal, DataInternal,...
                SigChargeDirections, [],[],Data{10},SigStructure,0};%, [], [], , Data{11}, , };
            
            U = SymTensor.SymMat2Ten(UMat,DataU);
            V = SymTensor.SymMat2Ten(VMat,DataV);
            Sig = SymTensor.SymMat2Ten(SigMat,DataSig);
            
            %{
            V.ChargeSide = -V.ChargeSide;
            V.ChargeDirections = -V.ChargeDirections;
            
            
            if ~FlagCorrectSides
                U.ChargeSide = -U.ChargeSide;
                U.ChargeDirections = -U.ChargeDirections;
                V.ChargeSide = -V.ChargeSide;
                V.ChargeDirections = -V.ChargeDirections;
                Sig.ChargeSide = -Sig.ChargeSide;
                Sig.ChargeDirections = -Sig.ChargeDirections;
            end
            
            %}
            
            
        end
        
        function self = Resize(self, NewChargeLegDimensions)
            %first modify ChargeLegDimension
            
            if ~isequal(size(self.ChargeLegDimensions),size(NewChargeLegDimensions))
                error('Error: The size of NewChargeLegDimension is incompatible with this Tensor')
            end
            
            if ~any(any(isnumeric(NewChargeLegDimensions) & IsInteger(NewChargeLegDimensions) & ~any(NewChargeLegDimensions<0)))
                error('Error: NewChargeLegDimension should only have positive integers (or zero) entries');
            end
            
            SizesOld = self.TensorEntriesSizes;
            SizesNew = zeros(size(SizesOld));
            SizesMid = zeros(size(SizesOld));
            for kk = 1:size(SizesNew,2)
                SizesNew(:,kk) = NewChargeLegDimensions(self.ChargeLabelsExternal(:,kk),kk);
                SizesMid(:,kk) = min(SizesNew(:,kk),SizesOld(:,kk));
            end
            self.ChargeLegDimensions = NewChargeLegDimensions;
            self.TensorEntriesSizes = SizesNew;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Cut out anything that is needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for kk = 1:numel(self.TensorEntries)
                NewSize = SizesMid(kk,:);
                OldSize = SizesOld(kk,:);
                if ~isequal(NewSize,OldSize)
                    Numbers = zeros(NewSize);
                    for aa = 1:size(NewSize,2)
                        Rep = NewSize;
                        Rep(aa) = 1;
                        Res = ones(size(NewSize));
                        Res(aa) = NewSize(aa);
                        Numbers = Numbers+repmat(reshape((0:(NewSize(aa)-1))*prod(OldSize(1:(aa-1))),Res),Rep);
                    end
                    Numbers = Numbers(:)+1;
                    
                    NewTensors = reshape(self.TensorEntries{kk}(Numbers),SizesMid(kk,:));
                    self.TensorEntries{kk} = NewTensors;
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %extend everything else with zeros
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for kk = 1:numel(self.TensorEntries)
                NewSize = SizesNew(kk,:);
                OldSize = SizesMid(kk,:);
                if ~isequal(NewSize,OldSize)
                    Numbers = zeros(OldSize);
                    for aa = 1:size(NewSize,2)
                        Rep = OldSize;
                        Rep(aa) = 1;
                        Res = ones(size(OldSize));
                        Res(aa) = OldSize(aa);
                        Numbers = Numbers+repmat(reshape((0:(OldSize(aa)-1))*prod(NewSize(1:(aa-1))),Res),Rep);
                    end
                    Numbers = Numbers(:)+1;
                    
                    NewTensors = zeros(SizesNew(kk,:));
                    NewTensors(Numbers) = self.TensorEntries{kk}(:);
                    self.TensorEntries{kk} = NewTensors;
                end
            end
            
            
        end
        
        function self = TenConj(self, Relabel,PerformConj)
            if nargin<2
                Relabel = false;
            end
            
            if nargin<3
                PerformConj = true;
            end
            
            
            if PerformConj
                for nn = 1:numel(self.TensorEntries)
                    self.TensorEntries{nn} = conj(self.TensorEntries{nn});
                end
            end
            
            self.ChargeDirections = -self.ChargeDirections;
            self.ChargeSide = -self.ChargeSide;
            self.ChargeSideInt = -self.ChargeSideInt;
            
            self.CyclicLegs = self.CyclicLegs(end:-1:1);
            self.CyclicLegsBraided = self.CyclicLegsBraided(end:-1:1);
            
            if Relabel
                self = self.BringToRegular(true);
            end
            
        end
        
        function self = TenConjNT(self, Relabel)
            if nargin<2
                Relabel = false;
            end
            
            for nn = 1:numel(self.TensorEntries)
                self.TensorEntries{nn} = conj(self.TensorEntries{nn});
            end
            
            %self.ChargeDirections = -self.ChargeDirections;
            %self.ChargeSide = -self.ChargeSide;
            %self.ChargeSideInt = -self.ChargeSideInt;
            
            %self.CyclicLegs = self.CyclicLegs(end:-1:1);
            %self.CyclicLegsBraided = self.CyclicLegsBraided(end:-1:1);
            
            if Relabel
                self = self.BringToRegular(true);
            end
        end
        
        function self = MakeHermitian(self)
            A = self.TenConj.BringToRegular(true);
            B = self.BringToRegular(true);
            
            if ~SymIsSameShape(A,B)
                error('Error: This tensor does not have a Hermitian shape when splitting from top to bottom (i.e. it is not a square matrix)');
            end
            self = (1/2)*(A+B);
            
        end
        
        function self = SortLabels(self)
            
            if ~self.SymHandle.IsNonAbelian
                %if it abelian
                AllDetailsRow = self.ChargeLabelsExternal;
                
            elseif ~self.SymHandle.IsMultiplicities
                %if non-abelian with no multiplicities
                AllDetailsRow = [self.ChargeLabelsExternal,self.ChargeLabelsInternal];
                
            else
                %if non-abelian with multiplicities
                AllDetailsRow = [self.ChargeLabelsExternal,self.ChargeLabelsInternal,self.MultiplicitiesInternal];
                
            end
            
            [AllDetailsRowMod,Index] = sortrows(AllDetailsRow);
            
            self.ChargeLabelsExternal = AllDetailsRowMod(:,1:size(self.ChargeLabelsExternal,2));
            
            if ~self.SymHandle.IsNonAbelian
                if ~isempty(self.ChargeLabelsInternal)
                    self.ChargeLabelsInternal = self.ChargeLabelsInternal(Index,:);
                end
            
            else
                if ~isempty(self.ChargeLabelsInternal)
                    self.ChargeLabelsInternal = AllDetailsRowMod(:,size(self.ChargeLabelsExternal,2)+(1:size(self.ChargeLabelsInternal,2)));
                end
                
                if self.SymHandle.IsMultiplicities
                    self.MultiplicitiesInternal = AllDetailsRowMod(:,size(self.ChargeLabelsExternal,2)+size(self.ChargeLabelsInternal,2)+(1:size(self.MultiplicitiesInternal,2)));
                    
                end
            end
            
            self.TensorEntries = self.TensorEntries(Index,:);
            self.TensorEntriesSizes = self.TensorEntriesSizes(Index,:);
            
        end
        
        function self = RandomSortLabels(self)
            Number = numel(self.TensorEntries);
            NewOrder = randperm(Number);
            
            self.ChargeLabelsExternal = self.ChargeLabelsExternal(NewOrder,:);
            self.ChargeLabelsInternal = self.ChargeLabelsInternal(NewOrder,:);
            self.MultiplicitiesInternal = self.MultiplicitiesInternal(NewOrder,:);
            self.TensorEntries = self.TensorEntries(NewOrder,:);
            self.TensorEntriesSizes = self.TensorEntriesSizes(NewOrder,:);
            
        end
        
        function self = SymReStructure(self, NewStructure, NewChargeDirections, NewChargeSideInt)
            self = self.SymReshape(NewStructure, NewChargeDirections, 1:numel(self.CyclicLegs), self.CapType, self.ChargeSide, NewChargeSideInt, ...
                self.Braidings, self.BraidingDirections, [],[],[],[]);
        end
        
        function self = SymReshape(self, NewStructure, NewChargeDirections, NewOrdering, NewCapType, NewChargeSide, NewChargeSideInt, ...
                NewBraidings, NewBraidingDirections, ModifyPhases, ExtraBraidings, ExtraBraidingDirections,RotateDirection)
            
            if nargin<13
                RotateDirection = 0;
            end
            %rotateDirection = +1  == clockwise
            %rotateDirection = -1  == counter clockwise
            %rotateDirection =  0  == Not set
            
            
            if self.FlagIsANumber
                %then there is no meaning to anything
                return;
            end
            
            %Inputs are:
            %NewStructure: obvious
            %NewChargeDirections: obvious
            %NewOrdering: what 
            %NewCapType: obvious
            %NewChargeSide: obvious
            %NewChargeSideInt: obvious
            %NewBraidings: obvious
            %NewBraidingsDirections: obvious
            %
            %ModifyPhases:n-by-2 matrix, first column is the label for the
            %leg in the original labeling, the second column tells us which
            %direction the modification is, the charge for pDown is
            %-MP(2)*Dir
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %checks:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isequal(sort(NewStructure(:)), sort(self.Structure(:)))
                error('Error: This structure is not compatible with the current structure')
            end
            
            if ~isequal(sort(NewOrdering(:)), (1:numel(self.ChargeSide))')
                error('Error:NewOrdering is not of the correct form')
            end
            
            
            DirChargesExternal = zeros([1,numel(self.ChargeSide)]);
            for kk = 1:numel(self.ChargeSide)
                DirChargesExternal(kk) = self.ChargeDirections(self.Structure == -kk);
            end
            
            
            DirChargesInternal = zeros([1,numel(self.ChargeSide)]);
            for kk = 1:numel(self.ChargeSideInt)
                DirChargesInternal(kk) = self.ChargeDirections(self.Structure == kk);
            end
            
            
            NewDirChargesExternal = zeros([1,numel(self.ChargeSide)]);
            for kk = 1:numel(NewChargeSide)
                NewDirChargesExternal(kk) = NewChargeDirections(NewStructure == -kk);
            end
            
            
            NewDirChargesInternal = zeros([1,numel(NewChargeSide)]);
            for kk = 1:numel(NewChargeSideInt)
                NewDirChargesInternal(kk) = NewChargeDirections(NewStructure == kk);
            end
            
            NewChargeSide = reshape(NewChargeSide,[1,numel(NewChargeSide)]);
            
            %Numbers = 1:numel(NewChargeSide);
            %CyclicLegsBraidedNew = Numbers(NewChargeSide == -1);
            %CyclicLegsBraidedNew = [Numbers(NewChargeSide == +1), CyclicLegsBraidedNew(end:-1:1)];
            
            %given just this data we need to start by matching the
            %legorders after the braiding takes place. Then our check is
            %that we must get -1... +1... for the
            %NewChargeSide(NewCyclicLegsBraided)
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%
            %Given all this Work out the new Cyclic legs
            %%%%%%%%%%%%%%%%%%%
            
            
            CyclicLegsOld = self.CyclicLegs;
            CyclicLegsBraidedOld = self.CyclicLegsBraided;
            OldBraidings = [self.Braidings(:);ExtraBraidings(:)];
            OldBraidingDirections = [self.BraidingDirections(:);ExtraBraidingDirections(:)];
            
            for kk = 1:numel(ExtraBraidings)
                gg = ExtraBraidings(kk);
                if gg~=numel(NewChargeSide)
                    CyclicLegsOld([gg+1,gg]) = CyclicLegsOld([gg,gg+1]);
                else
                    CyclicLegsOld([1,end]) = CyclicLegsOld([end,1]);
                end
            end
            
            %Now need to pull back to the innermost structure to work out
            %which one is up and which one is down
            
            CyclicLegsNew = CyclicLegsOld;
            
            CyclicLegsBraidedNew = CyclicLegsNew;
            for kk= 1:numel(NewBraidings)
                gg = NewBraidings(kk);
                if gg~=numel(NewChargeSide)
                    CyclicLegsBraidedNew([gg+1,gg]) = CyclicLegsBraidedNew([gg,gg+1]);
                else
                    CyclicLegsBraidedNew([1,end]) = CyclicLegsBraidedNew([end,1]);
                end
            end
            
            %this is the untransformed CyclicLegs
            
            %Now given this I can work out how much this needs to rotate to
            %agree with the structure:
            
            %find out which one is the first
            [~,Index] = sort(NewOrdering,'ascend');
            
            TopNumbersNew = NewChargeSide(Index) == +1;
            BottomNumbersNew = NewChargeSide(Index) == -1;
            TopNumbersOld = self.ChargeSide == +1;
            BottomNumbersOld = self.ChargeSide == -1;
            
            %-1 means top to bottom, +1 means bottom to top
            TopToBottomNames = (TopNumbersNew-BottomNumbersNew-TopNumbersOld+BottomNumbersOld)/2;
            
            %now in the order of CyclicLegsOld, CyclicLegs as we have
            %gotten rid of all the braiding by the point we impliment this part of the code.
            TopToBottom = TopToBottomNames(self.CyclicLegs); 
            
            
            NumberSwaps = sum(TopToBottom(1:(end-1)) ~= TopToBottom(2:end));
            
            if TopToBottom(1) == -1 && TopToBottom(end) == +1
                %edge case, either an error or 
                if NumberSwaps == 1
                    %if there are no zeros it doesn't matter which way this
                    %is done because of the lollypop identity
                    
                    AmountRotated = sum(TopToBottom == +1);
                elseif NumberSwaps == 2
                    
                    %if there are zeros then we can identify which side all
                    %the zeros are on (they should be the same)
                    
                    Value = self.ChargeSide(find(TopToBottomNames==0,1,'first'));
                    
                    if any(Value ~= self.ChargeSide(TopToBottomNames==0))
                        error('Error: The Rotations you have requested do not work');
                    end
                    
                    if Value == +1
                        %if the unbent ones are on the top then the first
                        %ones bend counterclockwise to the left
                        AmountRotated = -sum(TopToBottom == -1);
                    else %Value == -1
                        %if the unbent ones are on the bottom then the last
                        %ones bend clockwise to the left
                        AmountRotated = sum(TopToBottom == +1);
                    end
                else
                    error('Error: Too Many Swaps This Should Not Be Possible');
                end
                
                
            elseif TopToBottom(1) == -1
                %top to bottom, should be -ve amount rotated
                if NumberSwaps == 0
                    AmountRotated = 0;
                elseif NumberSwaps == 1
                    
                    if any(TopToBottomNames(self.ChargeSide == +1) ~=-1)
                        AmountRotated = -sum(TopToBottom == -1);
                    else
                        %special case
                        %this is difficult because all top legs are
                        %going to the bottom
                        
                        if RotateDirection == +1 %clockwise
                            AmountRotated = 0;
                        elseif RotateDirection == -1 %counterclockwise
                            AmountRotated = -sum(TopToBottom == -1);
                        else
                            error('Error: This SymReshape setting is ambiguous, please select clockwise or counter clockwise rotation. If this error doesn''t make sense to you please build the transformation out of other higher level SymTensor commands.')
                        end
                        
                        %Loc = find(CyclicLegsNew(1) == self.CyclicLegs, 1, 'first');
                        %AmountRotated = -(Loc);
                    end
                    
                elseif NumberSwaps < 4
                    Loc = find(TopToBottom ~= -1,1,'first');
                    if NumberSwaps == 3 && TopToBottom(Loc) == +1
                        error('Error: Wrong types of Swaps in SymReshape');
                    end
                    AmountRotated = -(Loc-1);
                else
                    error('Error: Too Many Swaps This Should Not Be Possible');
                end
                
                
                
            elseif TopToBottom(end) == +1
                %bottom to top, should be +ve amount rotated
                if NumberSwaps == 0
                    AmountRotated = 0;
                elseif NumberSwaps == 1
                    
                    if any (TopToBottomNames(self.ChargeSide == -1) ~=+1)
                        AmountRotated = sum(TopToBottom == +1);
                    else
                        %special case
                        %this is difficult because all bottom legs are
                        %going to the top
                        if RotateDirection == +1 %clockwise
                            AmountRotated = sum(TopToBottom == +1);
                        elseif RotateDirection == -1 %counterclockwise
                            AmountRotated = 0;
                        else
                            error('Error: This SymReshape setting is ambiguous, please select clockwise or counter clockwise rotation. If this error doesn''t make sense to you please build the transformation out of other higher level SymTensor commands.')
                        end
                        
                        %Loc = find(CyclicLegsNew(end) == self.CyclicLegs,1, 'first');
                        %AmountRotated = numel(TopToBottom)-Loc;
                    end
                    
                elseif NumberSwaps < 4
                    Loc = find(TopToBottom ~= +1,1,'last');
                    if NumberSwaps == 3 && TopToBottom(Loc) == -1
                        error('Error: Wrong types of Swaps in SymReshape');
                    end
                    AmountRotated = numel(TopToBottom)-Loc;
                    
                else
                    error('Error: Too Many Swaps This Should Not Be Possible');
                end
                
            else
                AmountRotated = 0;
            end

            
            %HERA  transform
            CyclicLegsBraidedNew = Index(CyclicLegsBraidedNew);
            CyclicLegsNew = Index(CyclicLegsNew);

                
            AmountRotatedPro = mod(AmountRotated, numel(CyclicLegsBraidedNew));
            
            CyclicLegsBraidedNew = CyclicLegsBraidedNew([(end-AmountRotatedPro+1):end,1:(end-AmountRotatedPro)]);
            CyclicLegsNew = CyclicLegsNew([(end-AmountRotatedPro+1):end,1:(end-AmountRotatedPro)]);
            
            Numbers = 1:numel(NewChargeSide);
            NewStoredLocations = cell([0,1]);
            for jj = max(NewStructure(:)):-1:1
                Updating = jj;
                Using = [];
                while ~isempty(Updating)
                    kk = Updating(end);
                    Updating(end) = [];
                    Use1 = NewStructure(1,kk+1);
                    Use2 = NewStructure(2,kk+1);
                    if Use1>0
                        if Use1>jj
                            Using = [Using,NewStoredLocations{Use1}];
                        else
                            Updating = [Updating,Use1];
                        end
                    else %Use1<0
                        Using = [Using, -Use1];
                    end
                    
                    if Use2>0
                        if Use2>jj
                            Using = [Using,NewStoredLocations{Use2}];
                        else
                            Updating = [Updating,Use2];
                        end
                    else %Use2<0
                        Using = [Using, -Use2];
                    end
                end
                NewStoredLocations{jj} = Using;
                
                
                %now checking, that these are all next to each other
                %note that Using has at least 2 elements
                
                WorkOn = Using(1);
                Using(1) = [];
                while ~isempty(Using)
                    CurrentUse = WorkOn(end);
                    
                    Temp = find(WorkOn(end) == CyclicLegsBraidedNew,1,'first');
                    WorkOn(end) = [];
                    
                    if numel(Temp)~=1
                        error('Affirmation Error: We should always be able to find this')
                    end
                    Temp = mod(Temp+[0,-2],Numbers(end))+1;
                    
                    TempFind1 = find(Using==CyclicLegsBraidedNew(Temp(1)));
                    
                    if numel(TempFind1)>1
                        error('Affirmation Error: We shouldn''t be able to find more then 1')
                    end
                    
                    if ~isempty(TempFind1)
                        WorkOn = [WorkOn,Using(TempFind1)];
                        Using(TempFind1) = [];
                    end
                    
                    TempFind2 = find(Using==CyclicLegsBraidedNew(Temp(2)));
                    
                    if numel(TempFind2)>1
                        error('Affirmation Error: We shouldn''t be able to find more then 1')
                    end
                    
                    if ~isempty(TempFind2)
                        WorkOn = [WorkOn,Using(TempFind2)];
                        Using(TempFind2) = [];
                    end
                    
                    if isempty(WorkOn)
                        error('Error: This is not a permissable tensor');
                    end
                end
                
            end
            
            
            
            
            
            
            
            
            
            
            %HERA  transform
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{kk} = permute(self.TensorEntries{kk},NewOrdering);
            end
            
            
            
            
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STEP 0) work out Caps
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STEP 1) work out braidings for new and old
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ChargeDirectionsExternalOld = DirChargesExternal(CyclicLegsBraidedOld);
            ChargeDirectionsExternalNew = NewDirChargesExternal(CyclicLegsBraidedNew);
            
            MonoBraiding = [];
            TwistType = [];
            ExtraRounds = []; %+1 means right to left(Round), and -1 means left to right (RoundInv)
            
            %%%%%%%%%%%%%%%%%%%
            %OldBraidings
            %%%%%%%%%%%%%%%%%%%
            for kk = 1:numel(OldBraidings)
                ll = OldBraidings(kk);
                if ll~=numel(self.ChargeSide)
                    MonoBraiding = [MonoBraiding, ll];
                    TwistType = [TwistType,OldBraidingDirections(kk)];
                    ExtraRounds = [ExtraRounds, 0];
                elseif OldBraidingDirections(kk)==+1
                    MonoBraiding = [MonoBraiding, (numel(self.ChargeSide)-1):-1:2,1:(numel(self.ChargeSide)-1)];
                    TwistType = [TwistType,-ones([1,numel(self.ChargeSide)-2]),ones([1,numel(self.ChargeSide)-1])];
                    if numel(self.ChargeSide)>2
                        ExtraRounds = [ExtraRounds,+1,zeros([1,numel(self.ChargeSide)-3]),-1,zeros([1,numel(self.ChargeSide)-2])];
                    else
                        ExtraRounds = [ExtraRounds,1,0];
                    end
                else %OldBraidingDirections(kk) == -1
                    MonoBraiding = [MonoBraiding, 1:(numel(self.ChargeSide)-2),(numel(self.ChargeSide)-1):-1:1];
                    TwistType = [TwistType,ones([1,numel(self.ChargeSide)-2]),-ones([1,numel(self.ChargeSide)-1])];
                    if numel(self.ChargeSide)>2
                        ExtraRounds = [ExtraRounds,-1,zeros([1,numel(self.ChargeSide)-3]),1,zeros([1,numel(self.ChargeSide)-2])];
                    else
                        ExtraRounds = [ExtraRounds,1,0];
                    end
                end
            end
            
            BraidingOld = cell([0,2]);
            for kk = 1:numel(MonoBraiding)
                gg = MonoBraiding(kk);
                if ExtraRounds(kk) == -1
                    BraidingOld = [BraidingOld;{'RoundInv',[-numel(self.ChargeSide);ChargeDirectionsExternalOld(numel(self.ChargeSide))]}];
                elseif ExtraRounds(kk) == +1
                    BraidingOld = [BraidingOld;{'Round',[-1;ChargeDirectionsExternalOld(1)]}];
                end
                cLabel = -gg;
                cDir = ChargeDirectionsExternalOld(gg);
                
                bLabel = -(gg+1);
                bDir = ChargeDirectionsExternalOld(gg+1);
                
                if gg == numel(ChargeDirectionsExternalOld)-1
                    dLabel = nan;
                    dDir = 0;
                else
                    dLabel = gg;
                    dDir = +1;
                end
                
                eLabel = gg-1;
                eDir = +1;
                MultLabel1 = eLabel;
                MultLabel2 = dLabel;
                
                %this will only occur if it is a B type
                if gg == 2
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalOld(1);
                else
                    aLabel = gg-2;
                    aDir = +1;
                end
                
                if TwistType(kk) == 1 %(Left)
                    if gg == 1
                        BraidingOld = [BraidingOld; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                    else
                        BraidingOld = [BraidingOld; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                    end
                else %TwistType(kk) == -1 %(Right)
                    if gg == 1
                        BraidingOld = [BraidingOld; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                    else
                        BraidingOld = [BraidingOld; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                    end
                end
                ChargeDirectionsExternalOld([gg+1,gg]) = ChargeDirectionsExternalOld([gg,gg+1]);
            end
            
            MonoBraiding = [];
            TwistType = [];
            ExtraRounds = []; %+1 means right to left(Round), and -1 means left to right (RoundInv)
            
            
            %%%%%%%%%%%%%%%%%%%
            %NewBraidings
            %%%%%%%%%%%%%%%%%%%
            for kk = 1:numel(NewBraidingDirections)
                ll = NewBraidings(kk);
                if ll~=numel(NewChargeSide)
                    MonoBraiding = [MonoBraiding, ll];
                    TwistType = [TwistType,NewBraidingDirection(kk)];
                    ExtraRounds = [ExtraRounds, 0];
                elseif NewBraidingDirection(kk)==+1
                    MonoBraiding = [MonoBraiding, (numel(NewChargeSide)-1):-1:2,1:(numel(NewChargeSide)-1)];
                    TwistType = [TwistType,-ones([1,numel(NewChargeSide)-2]),ones([1,numel(NewChargeSide)-1])];
                    if numel(self.ChargeSide)>2
                        ExtraRounds = [ExtraRounds,+1,zeros([1,numel(NewChargeSide)-3]),-1,zeros([1,numel(NewChargeSide)-2])];
                    else
                        ExtraRounds = [ExtraRounds,1,0];
                    end
                else %NewBraidingDirection(kk) == -1
                    MonoBraiding = [MonoBraiding, 1:(numel(NewChargeSide)-2),(numel(NewChargeSide)-1):-1:1];
                    TwistType = [TwistType,ones([1,numel(NewChargeSide)-2]),-ones([1,numel(NewChargeSide)-1])];
                    if numel(self.ChargeSide)>2
                        ExtraRounds = [ExtraRounds,-1,zeros([1,numel(NewChargeSide)-3]),1,zeros([1,numel(NewChargeSide)-2])];
                    else
                        ExtraRounds = [ExtraRounds,1,0];
                    end
                end
            end
            
            BraidingNew = cell([0,2]);
            for kk = 1:numel(MonoBraiding)
                gg = MonoBraiding(kk);
                if ExtraRounds(kk) == -1
                    BraidingNew = [BraidingNew;{'RoundInv',[-numel(self.ChargeSide);ChargeDirectionsExternalNew(numel(self.ChargeSide))]}];
                elseif ExtraRounds(kk) == +1
                    BraidingNew = [BraidingNew;{'Round',[-1;ChargeDirectionsExternalNew(1)]}];
                end
                cLabel = -gg;
                cDir = ChargeDirectionsExternalNew(gg);
                
                bLabel = -(gg+1);
                bDir = ChargeDirectionsExternalNew(gg+1);
                
                if gg == numel(ChargeDirectionsExternalOld)-1
                    dLabel = nan;
                    dDir = 0;
                else
                    dLabel = gg;
                    dDir = +1;
                end
                
                eLabel = gg-1;
                eDir = +1;
                MultLabel1 = eLabel;
                MultLabel2 = dLabel;
                
                %this will only occur if it is a B type
                if gg == 2
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalNew(1);
                else
                    aLabel = gg-2;
                    aDir = +1;
                end
                
                if TwistType(kk) == 1 %(Left)
                    if gg == 1
                        BraidingNew = [BraidingNew; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                    else
                        BraidingNew = [BraidingNew; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                    end
                else %TwistType(kk) == -1 %(Right)
                    if gg == 1
                        BraidingNew = [BraidingNew; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                    else
                        BraidingNew = [BraidingNew; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                    end
                end
                ChargeDirectionsExternalNew([gg+1,gg]) = ChargeDirectionsExternalNew([gg,gg+1]);
            end

            ReverseBraidings = BraidingNew(end:-1:1,:);
            GenerateBraidings = BraidingOld;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STEP 2)  work out phases
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %ignore for the moment
            GeneratePhases = cell([0,2]);
            ReversePhases = cell([0,2]);
            
            for kk = 1:size(ModifyPhases,1)
                gg = find(CyclicLegsBraidedOld == ModifyPhases(kk,1));
                GeneratePhases = [GeneratePhases;repmat({'pDown', [-gg;sign(ModifyPhases(kk,2))*DirChargesExternal(gg)]},[abs(ModifyPhases(kk,2)),1])];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STEP 3) work out rotate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ConvertRotations = cell([0,2]);
            MonoBraiding = [];
            TwistType = [];
            if AmountRotated<0
                for kk = 1:abs(AmountRotated)
                    MonoBraiding = [MonoBraiding, 1:(numel(NewChargeSide)-1)];
                    TwistType = [TwistType, +ones([1,numel(NewChargeSide)-1])];
                    %the first term makes the bottom an alternative, the second
                    %does the rounding.
                    ConvertRotations = [ConvertRotations; {'pDown',[-kk;-ChargeDirectionsExternalOld(kk)]};{'RoundInv', [-kk;ChargeDirectionsExternalOld(kk)]}];
                end
            else %AmountRotated>0 % this means last goes to first location
                for kk = numel(NewChargeSide):-1:(numel(NewChargeSide)-AmountRotated+1)
                    MonoBraiding = [MonoBraiding, (numel(NewChargeSide)-1):-1:1];
                    TwistType = [TwistType, -ones([1,numel(NewChargeSide)-1])];
                    %the first term makes the bottom an alternative, the second
                    %does the rounding.
                    ConvertRotations = [ConvertRotations; {'pDown',[-kk;ChargeDirectionsExternalOld(kk)]};{'pDown',[-kk;ChargeDirectionsExternalOld(kk)]};{'Round', [-kk;ChargeDirectionsExternalOld(kk)]}];
                end
            end
            
            for kk = 1:numel(MonoBraiding)
                gg = MonoBraiding(kk);
                cLabel = -gg;
                cDir = ChargeDirectionsExternalOld(gg);
                
                bLabel = -(gg+1);
                bDir = ChargeDirectionsExternalOld(gg+1);
                
                if gg == numel(ChargeDirectionsExternalOld)-1
                    dLabel = nan;
                    dDir = 0;
                else
                    dLabel = gg;
                    dDir = +1;
                end
                
                eLabel = gg-1;
                eDir = +1;
                MultLabel1 = eLabel;
                MultLabel2 = dLabel;
                
                %this will only occur if it is a B type
                if gg == 2
                    aLabel = -1;
                    aDir = ChargeDirectionsExternalOld(1);
                else
                    aLabel = gg-2;
                    aDir = +1;
                end
                
                if TwistType(kk) == 1 %(Left)
                    if gg == 1
                        ConvertRotations = [ConvertRotations; {'RInv',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                    else
                        ConvertRotations = [ConvertRotations; {'BInv',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                    end
                else %TwistType(kk) == -1 %(Right)
                    if gg == 1
                        ConvertRotations = [ConvertRotations; {'R',[bLabel,cLabel,dLabel,MultLabel2,MultLabel2;bDir,cDir,dDir,0,0]}];%disp('RInv0')
                    else
                        ConvertRotations = [ConvertRotations; {'B',[aLabel,bLabel,cLabel,dLabel,eLabel,eLabel,MultLabel1,MultLabel2,MultLabel1,MultLabel2;aDir,bDir,cDir,dDir,eDir,eDir,0,0,0,0]}]; %disp('BInv1')
                    end
                end
                ChargeDirectionsExternalOld([gg+1,gg]) = ChargeDirectionsExternalOld([gg,gg+1]);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STEP 4) use from standard
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            GenerateFromStandard = SymTensor.FromStandard(self,1);
            %rotations have been taken care of by the ConvertRotations
            %section.
            ReverseFromStandard = SymTensor.FromStandard({NewStructure, NewChargeDirections, NewChargeSide, NewChargeSideInt, CyclicLegsBraidedNew},1,true);
            
            NumberExternal = numel(CyclicLegsBraidedNew);
            NumberInternal = max(0,NumberExternal-2);
            [~,ReverseCyclicLegsBraidedNew] = sort(CyclicLegsBraidedNew,'ascend');
            
            FixSwapReverse = [-CyclicLegsBraidedNew,1:NumberInternal;...
                zeros(1,NumberExternal+NumberInternal);...
                1:NumberInternal,zeros(1,NumberExternal);...
                -ReverseCyclicLegsBraidedNew,1:NumberInternal;...
                zeros(1,NumberExternal+NumberInternal);...
                1:NumberInternal,zeros(1,NumberExternal)];
            
            [~,ReverseCyclicLegsBraidedOld] = sort(CyclicLegsBraidedOld,'ascend');
            
            %note that this reverses the swap that FromStandard does, the
            %other one had to impliment the stanard of going from -1 to -N
            FixSwap = [-ReverseCyclicLegsBraidedOld,1:NumberInternal;...
                zeros(1,NumberExternal+NumberInternal);...
                1:NumberInternal,zeros(1,NumberExternal);...
                -CyclicLegsBraidedOld,1:NumberInternal;...
                zeros(1,NumberExternal+NumberInternal);...
                1:NumberInternal,zeros(1,NumberExternal)];
            
            
            if ~isempty(ReverseBraidings)
                for kk = 1:size(ReverseBraidings,2)
                    if isequal(ReverseBraidings{kk,1}, 'R')
                        ReverseBraidings{kk,1} = 'RInv';
                    elseif isequal(ReverseBraidings{kk,1}, 'RInv')
                        ReverseBraidings{kk,1} = 'R';
                    elseif isequal(ReverseBraidings{kk,1}, 'B')
                        ReverseBraidings{kk,1} = 'BInv';
                    elseif isequal(ReverseBraidings{kk,1}, 'BInv')
                        ReverseBraidings{kk,1} = 'B';
                    elseif isequal(ReverseBraidings{kk,1}, 'Round')
                        ReverseBraidings{kk,1} = 'RoundInv';
                    elseif isequal(ReverseBraidings{kk,1}, 'RoundInv')
                        ReverseBraidings{kk,1} = 'Round';
                    else
                        error('Affirmation Error: there is something other then R, RInv, B, BInv, Round, RoundInv in ReverseBraidings')
                    end
                end
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %STEP 5) Now convert
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %add together
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            FullGenerateWords = [GenerateFromStandard; GeneratePhases; GenerateBraidings;...
                ConvertRotations; ReverseBraidings; ReversePhases; ReverseFromStandard];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %generate numbers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if self.SymHandle.IsNonAbelian
                ConvertNumbers = self.SymHandle.GenerateNumbers(FullGenerateWords,...
                    self.ChargeLabelsExternal,self.ChargeLabelsInternal,self.MultiplicitiesInternal);
            else
                ConvertNumbers = self.SymHandle.GenerateNumbers(FullGenerateWords,...
                    self.ChargeLabelsExternal,self.ChargeLabelsInternal,self.MultiplicitiesInternal);
%                ConvertNumbers = self.SymHandle.GenerateNumbersAbelian(FullGenerateWords,self.ChargeLabelsExternal);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %transform
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            LabelsIn = ConvertNumbers{1};
            LabelsOut = ConvertNumbers{2};
            [ListIn,ListOut,Entries] = find(ConvertNumbers{3});
            
            [~,A,C] = unique([LabelsIn;[self.ChargeLabelsExternal,self.ChargeLabelsInternal,self.MultiplicitiesInternal]], 'rows','first');
            %[~,Index] = sort(A,'ascend');
            Index = A(C);
            Index = Index((size(LabelsIn,1)+1):end);
            ListIn = Index(ListIn);
            
            NewChargeLabelsExternal = LabelsOut(:,1:size(self.ChargeLabelsExternal,2));
            NewChargeLabelsInternal = LabelsOut(:,size(NewChargeLabelsExternal,2) + (1:size(self.ChargeLabelsInternal,2)));
            NewMultiplicitiesInternal = LabelsOut(:,size(NewChargeLabelsExternal,2) + size(NewChargeLabelsInternal,2) + (1:size(self.MultiplicitiesInternal,2)));
            
            NewTensorEntries = cell([max(ListOut),1]);
            NewTensorEntriesSizes = zeros([max(ListOut),size(self.TensorEntriesSizes,2)]);
            FilledEntry = false(size(NewTensorEntries));
            for kk = 1:numel(ListIn)
                if FilledEntry(ListOut(kk))
                    NewTensorEntries{ListOut(kk)} = NewTensorEntries{ListOut(kk)} + Entries(kk)*self.TensorEntries{ListIn(kk)};
                else
                    FilledEntry(ListOut(kk)) = true;
                    NewTensorEntries{ListOut(kk)} = Entries(kk)*self.TensorEntries{ListIn(kk)};
                    NewTensorEntriesSizes(ListOut(kk),:) = self.TensorEntriesSizes(ListIn(kk),:);
                end
            end
            
            NewTensorEntries = NewTensorEntries(FilledEntry);
            NewTensorEntriesSizes = NewTensorEntriesSizes(FilledEntry,NewOrdering);
            self.ChargeLegDimensions = self.ChargeLegDimensions(:,NewOrdering);
            
            %test HERA
            self.ChargeLabelsExternal = NewChargeLabelsExternal(FilledEntry,:);%NewOrdering);
            self.ChargeLabelsInternal = NewChargeLabelsInternal(FilledEntry,:);
            self.MultiplicitiesInternal = NewMultiplicitiesInternal(FilledEntry,:);
            
            self.TensorEntries = NewTensorEntries;
            self.TensorEntriesSizes = NewTensorEntriesSizes;
            
            self.Structure = NewStructure;
            self.ChargeDirections = NewChargeDirections;
            
            self.Braidings = NewBraidings;
            self.BraidingDirections = NewBraidingDirections;
            self.CapType = NewCapType;
            
            self.ChargeSide = NewChargeSide;
            self.ChargeSideInt = NewChargeSideInt;
            
            self.CyclicLegs = CyclicLegsNew;
            self.CyclicLegsBraided = CyclicLegsBraidedNew;
            
            self.StoredLocations = NewStoredLocations;
            
        end
        
        function self = IdentityToRight(self, ChargesDims)
            
            DimNumber = size(self.SymHandle.Dim,1);
            TrivialCharges = self.SymHandle.TrivialIrrep;
            
            if ~isnumeric(ChargesDims)
                error('Error: ChargesDim should be a numeric number')
            end
            
            if numel(ChargesDims) == 1
                ChargesDims = repmat(ChargesDims,[DimNumber,1]);
            elseif numel(ChargesDims) == DimNumber
                ChargesDims = reshape(ChargesDims,[DimNumber,1]);
            elseif numel(ChargesDims) < DimNumber
                ChargesDims = [reshape(ChargesDims,[numel(ChargesDims),1]);zeros([DimNumber-numel(ChargesDims),1])];
            else
                error('Error: For this symmetry there should be less a different number of DimNumbers');
            end
            
            if ~IsInteger(ChargesDims) || any(ChargesDims+10^-14)<0
                error('Error: ChargesDim should be positive semi-definite integers only')
            end
            ChargesDims = round(ChargesDims);
            DimNumbersGood = sum(ChargesDims>0);
            
            if DimNumbersGood == 0
                error('Error: we should have at least one ChargesDim which is non-zero');
            end
            
            if self.FlagIsANumber
                self = SymTensor.CreateIdentity(self.SymHandle,ChargesDims,ChargesDims, [1,1])*self.getNumber;
                return;
            end
            TopNumber = sum(self.ChargeSide == +1);
            BottomNumber = sum(self.ChargeSide == -1);
            AllNumber = TopNumber+BottomNumber;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fix braidings
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            self.Braidings = self.Braidings+2*(self.Braidings>TopNumber);
            if any(self.Braidings == TopNumber)
                error('So far I can''t add an identity over a braid')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fix CapType and CyclicLegs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %need to adjust for caps and cups
            self.CapType = [self.CapType(1:TopNumber),0,0,self.CapType((TopNumber+1):end)];
            self.CyclicLegs = [self.CyclicLegs(1:TopNumber),AllNumber+[2,1],self.CyclicLegs((TopNumber+1):end)];
            self.CyclicLegsBraided = [self.CyclicLegsBraided(1:TopNumber),AllNumber+[2,1],self.CyclicLegsBraided((TopNumber+1):end)];
            %self.FreePermuteTracker = [self.FreePermuteTracker,AllNumber+[1,2]];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fix structure and ChargeDirections (and ChargeSide)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            self.ChargeLegDimensions = [self.ChargeLegDimensions,ChargesDims,ChargesDims];
            
            if TopNumber>0
                CoreConnected = -self.CyclicLegsBraided(TopNumber);
                
                self.ChargeDirections = [self.ChargeDirections,[self.ChargeDirections(self.Structure == CoreConnected),+1;+1,-1]];
                
                self.Structure(self.Structure == CoreConnected) = AllNumber-1;
                self.Structure = [self.Structure,[CoreConnected, -(AllNumber+2); AllNumber, -(AllNumber+1)]];
                
                self.ChargeSideInt = [self.ChargeSideInt, +1,+1];
                self.ChargeSide = [self.ChargeSide, -1,+1];
                
            else
                %then has to be done at the bottom
                CoreConnected = -self.CyclicLegsBraided(TopNumber+3); %+3 because we have to take into account the added legs
                
                self.ChargeDirections = [self.ChargeDirections,[self.ChargeDirections(self.Structure == CoreConnected),+1;+1,-1]];
                
                self.Structure(self.Structure == CoreConnected) = AllNumber-1;
                self.Structure = [self.Structure,[CoreConnected, -(AllNumber+2); AllNumber, -(AllNumber+1)]];
                
                self.ChargeSideInt = [self.ChargeSideInt, +1,-1];
                self.ChargeSide = [self.ChargeSide, -1,+1];
                
            end
            
            for kk = 1:numel(self.StoredLocations)
                if any(self.StoredLocations{kk} == -CoreConnected)
                    self.StoredLocations{kk} = [self.StoredLocations{kk},AllNumber+[1,2]];
                end
            end
            
            self.StoredLocations{AllNumber-2+1} = [-CoreConnected, AllNumber+[1,2]];
            self.StoredLocations{AllNumber-2+2} = AllNumber+[1,2];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %update ChargeLabels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ~isempty(self.MultiplicitiesInternal)
                self.MultiplicitiesInternal = [repmat(self.MultiplicitiesInternal,[DimNumbersGood,1]),ones([size(self.MultiplicitiesInternal,1)*DimNumbersGood,2])];
                self.ChargeLabelsInternal = repmat([self.ChargeLabelsInternal, self.ChargeLabelsExternal(:,-CoreConnected),...
                                        repmat(TrivialCharges,[size(self.ChargeLabelsExternal,1),1])],[DimNumbersGood,1]);
            else
                self.MultiplicitiesInternal = ones([size(self.ChargeLabelsExternal,1)*DimNumbersGood,2]);
                self.ChargeLabelsInternal = [repmat(self.ChargeLabelsExternal(:,-CoreConnected),[DimNumbersGood,1]),...
                                        repmat(TrivialCharges,[size(self.ChargeLabelsExternal,1)*DimNumbersGood,1])];
            end
            
            Numbers = 1:DimNumber;
            self.ChargeLabelsExternal = [repmat(self.ChargeLabelsExternal,[DimNumbersGood,1]),...
                reshape(repmat(reshape(Numbers(ChargesDims ~=0),[1,1,DimNumbersGood]),[2,size(self.ChargeLabelsExternal,1)]),[2,DimNumbersGood*size(self.ChargeLabelsExternal,1)])'];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update TensorEntries
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            self.TensorEntries = repmat(self.TensorEntries,[DimNumbersGood,1]);
            
            GoodSizesNew = ChargesDims(ChargesDims ~= 0);
            self.TensorEntriesSizes = [repmat(self.TensorEntriesSizes,[DimNumbersGood,1]),...
                reshape(repmat(reshape(GoodSizesNew,[1,1,DimNumbersGood]),[2,size(self.TensorEntriesSizes,1)]),[2,DimNumbersGood*size(self.TensorEntriesSizes,1)])'];
            
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{kk} = reshape(kron(reshape(self.TensorEntries{kk},[numel(self.TensorEntries{kk}),1]),...
                    reshape(eye(self.TensorEntriesSizes(kk,end)),[1,self.TensorEntriesSizes(kk,end)])),self.TensorEntriesSizes(kk,:));
            end
            
        end
        
        function self = IdentityToLeft(self, ChargesDims)
            
            DimNumber = size(self.SymHandle.Dim,1);
            TrivialCharges = self.SymHandle.TrivialIrrep;
            
            if ~isnumeric(ChargesDims)
                error('Error: ChargesDim should be a numeric number')
            end
            
            if numel(ChargesDims) == 1
                ChargesDims = repmat(ChargesDims,[DimNumber,1]);
            elseif numel(ChargesDims) == DimNumber
                ChargesDims = reshape(ChargesDims,[DimNumber,1]);
            elseif numel(ChargesDims) < DimNumber
                ChargesDims = [reshape(ChargesDims,[numel(ChargesDims),1]);zeros([DimNumber-numel(ChargesDims),1])];
            else
                error('Error: For this symmetry there should be less a different number of DimNumbers');
            end
            
            if ~IsInteger(ChargesDims) || any(ChargesDims+10^-14)<0
                error('Error: ChargesDim should be positive semi-definite integers only')
            end
            ChargesDims = round(ChargesDims);
            DimNumbersGood = sum(ChargesDims>0);
            
            if DimNumbersGood == 0
                error('Error: we should have at least one ChargesDim which is non-zero');
            end
            
            if self.FlagIsANumber
                self = SymTensor.CreateIdentity(self.SymHandle,ChargesDims,ChargesDims, [1,1])*self.getNumber;
                return;
            end
            TopNumber = sum(self.ChargeSide == +1);
            BottomNumber = sum(self.ChargeSide == -1);
            AllNumber = TopNumber+BottomNumber;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fix braidings
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            self.Braidings = self.Braidings+1+1*(self.Braidings==AllNumber);
            if any(self.Braidings == AllNumber)
                error('So far I can''t add an identity over a braid')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fix CapType and CyclicLegs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %need to adjust for caps and cups
            self.CapType = [self.CapType(1:TopNumber),0,0,self.CapType((TopNumber+1):end)];
            self.CyclicLegs = [AllNumber+2,self.CyclicLegs,AllNumber+1];
            self.CyclicLegsBraided = [AllNumber+2,self.CyclicLegsBraided,AllNumber+1];
            %self.FreePermuteTracker = [self.FreePermuteTracker,AllNumber+[1,2]];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %fix structure and ChargeDirections (and ChargeSide)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            self.ChargeLegDimensions = [self.ChargeLegDimensions,ChargesDims,ChargesDims];
            
            if TopNumber>0
                CoreConnected = -self.CyclicLegsBraided(2); %2 because of the identity
                
                self.ChargeDirections = [self.ChargeDirections,[self.ChargeDirections(self.Structure == CoreConnected),+1;+1,-1]];
                
                self.Structure(self.Structure == CoreConnected) = AllNumber-1;
                self.Structure = [self.Structure,[CoreConnected, -(AllNumber+2); AllNumber, -(AllNumber+1)]];
                
                self.ChargeSideInt = [self.ChargeSideInt, +1,+1];
                self.ChargeSide = [self.ChargeSide, -1,+1];
                
            else
                %then has to be done at the bottom
                CoreConnected = -self.CyclicLegsBraided(end-1); %-1 because we have to take into account the added legs
                
                self.ChargeDirections = [self.ChargeDirections,[self.ChargeDirections(self.Structure == CoreConnected),+1;+1,-1]];
                
                self.Structure(self.Structure == CoreConnected) = AllNumber-1;
                self.Structure = [self.Structure,[CoreConnected, -(AllNumber+2); AllNumber, -(AllNumber+1)]];
                
                self.ChargeSideInt = [self.ChargeSideInt, +1,-1];
                self.ChargeSide = [self.ChargeSide, -1,+1];
                
            end
            
            for kk = 1:numel(self.StoredLocations)
                if any(self.StoredLocations{kk} == -CoreConnected)
                    self.StoredLocations{kk} = [self.StoredLocations{kk},AllNumber+[1,2]];
                end
            end
            
            self.StoredLocations{AllNumber-2+1} = [-CoreConnected, AllNumber+[1,2]];
            self.StoredLocations{AllNumber-2+2} = AllNumber+[1,2];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %update ChargeLabels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(self.MultiplicitiesInternal)
                self.MultiplicitiesInternal = [repmat(self.MultiplicitiesInternal,[DimNumbersGood,1]),ones([size(self.MultiplicitiesInternal,1)*DimNumbersGood,2])];
                self.ChargeLabelsInternal = repmat([self.ChargeLabelsInternal, self.ChargeLabelsExternal(:,-CoreConnected),...
                                        repmat(TrivialCharges,[size(self.ChargeLabelsExternal,1),1])],[DimNumbersGood,1]);
            else
                self.MultiplicitiesInternal = ones([size(self.ChargeLabelsExternal,1)*DimNumbersGood,2]);
                self.ChargeLabelsInternal = [repmat(self.ChargeLabelsExternal(:,-CoreConnected),[DimNumbersGood,1]),...
                                        repmat(TrivialCharges,[size(self.ChargeLabelsExternal,1)*DimNumbersGood,1])];
            end
            
            
            Numbers = 1:DimNumber;
            self.ChargeLabelsExternal = [repmat(self.ChargeLabelsExternal,[DimNumbersGood,1]),...
                reshape(repmat(reshape(Numbers(ChargesDims ~=0),[1,1,DimNumbersGood]),[2,size(self.ChargeLabelsExternal,1)]),[2,DimNumbersGood*size(self.ChargeLabelsExternal,1)])'];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update TensorEntries
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            self.TensorEntries = repmat(self.TensorEntries,[DimNumbersGood,1]);
            
            GoodSizesNew = ChargesDims(ChargesDims ~= 0);
            self.TensorEntriesSizes = [repmat(self.TensorEntriesSizes,[DimNumbersGood,1]),...
                reshape(repmat(reshape(GoodSizesNew,[1,1,DimNumbersGood]),[2,size(self.TensorEntriesSizes,1)]),[2,DimNumbersGood*size(self.TensorEntriesSizes,1)])'];
            
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{kk} = reshape(kron(reshape(self.TensorEntries{kk},[numel(self.TensorEntries{kk}),1]),...
                    reshape(eye(self.TensorEntriesSizes(kk,end)),[1,self.TensorEntriesSizes(kk,end)])),self.TensorEntriesSizes(kk,:));
            end
            
        end
        
        function self = ReOrderLegs(self, NewOrder)
            if ~isequal(1:numel(self.ChargeSide),sort(NewOrder(:)','ascend'))
                error('Error: NewOrder should have the same number of entries as the number of legs and count from 1 to the maximum')
            end
            
            [~,Index] = sort(NewOrder,'ascend');
            %Temp = self.ChargeDirections(self.Structure<0);
            %self.ChargeDirections(self.Structure<0) = Temp(NewOrder);
            self.Structure(self.Structure<0) = -Index(-self.Structure(self.Structure<0));
            
            
            
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{kk} = permute(self.TensorEntries{kk},NewOrder);
            end
            self.TensorEntriesSizes = self.TensorEntriesSizes(:,NewOrder);
            self.ChargeLabelsExternal = self.ChargeLabelsExternal(:,NewOrder);
            for kk = 1:numel(self.StoredLocations)
                self.StoredLocations{kk} = NewOrder(self.StoredLocations{kk});
            end
            self.CapType = self.CapType(NewOrder);
            self.ChargeLegDimensions = self.ChargeLegDimensions(:,NewOrder);
            self.ChargeSide = self.ChargeSide(NewOrder);
            self.CyclicLegs = Index(self.CyclicLegs);
            self.CyclicLegsBraided = Index(self.CyclicLegsBraided);
            
        end
        
        function self = FlipLegDirections(self, FlippedLegs)
            for kk = 1:numel(FlippedLegs)
                self.ChargeDirections(FlippedLegs(kk)==self.Structure) = -self.ChargeDirections(FlippedLegs(kk)==self.Structure);
            end
            self.ChargeLabelsExternal(:,-FlippedLegs(FlippedLegs(:)<0)) = self.SymHandle.InverseIrrep(self.ChargeLabelsExternal(:,-FlippedLegs(FlippedLegs(:)<0)));
            self.ChargeLabelsInternal(:,FlippedLegs(FlippedLegs(:)>0)) = self.SymHandle.InverseIrrep(self.ChargeLabelsInternal(:,FlippedLegs(FlippedLegs(:)>0)));
        end
        
        function self = BringToRegular(self,RenameLegs)
            
            if nargin<2
                RenameLegs = false;
            end
            
            %this brings it to a point where there are no braids and it is
            %an hourglass structure just and would be created by
            %SymTensor.CreateRandomUnitary
            
            Numbers = 1:numel(self.ChargeSide);
            TopLegsOriginal = Numbers(self.ChargeSide == +1)';
            TopLegsOriginal = self.CyclicLegs(any(repmat(self.CyclicLegs, [numel(TopLegsOriginal),1]) == repmat(TopLegsOriginal, [1,numel(self.CyclicLegs)]),1));
            
            BottomLegsOriginal = Numbers(self.ChargeSide == -1)';
            BottomLegsOriginal = self.CyclicLegs(any(repmat(self.CyclicLegs, [numel(BottomLegsOriginal),1]) == repmat(BottomLegsOriginal, [1,numel(self.CyclicLegs)]),1));
            
            CyclicLegsOriginal = [TopLegsOriginal, BottomLegsOriginal];
            %this is the order that I need the new legs to be in
            
            BottomLegsOriginal = BottomLegsOriginal(end:-1:1);
            
            ExternalChargeDirections = zeros([numel(CyclicLegsOriginal),1]);
            for kk = 1:numel(CyclicLegsOriginal)
                ExternalChargeDirections(kk) = self.ChargeDirections(self.Structure == -kk);
            end
            
            TotalLegsBottom = numel(BottomLegsOriginal);
            TotalLegsTop = numel(TopLegsOriginal);
            ChargeDirectionsUsedBottom = ExternalChargeDirections(BottomLegsOriginal);
            ChargeDirectionsUsedTop = ExternalChargeDirections(TopLegsOriginal);
            
            if numel(CyclicLegsOriginal)>=2
               NewStructure = zeros([2,numel(CyclicLegsOriginal)-1]);
               NewStructure = ones([2,numel(CyclicLegsOriginal)-1]);
               if TotalLegsBottom*TotalLegsTop >0
                    NewStructure = zeros([2,TotalLegsBottom+TotalLegsTop-1]);
                    NewChargeDirections = ones([2,TotalLegsBottom+TotalLegsTop-1]);
                    
                    if TotalLegsBottom == 1
                        NewStructure(1,1) = -BottomLegsOriginal(1);
                        NewChargeDirections(1,1) = ChargeDirectionsUsedBottom(1);
                    else
                        NewStructure(1,1) = TotalLegsBottom-1;
                        NewChargeDirections(1,1) = -1;
                        NewStructure(2,TotalLegsBottom) = -BottomLegsOriginal(end);
                        NewChargeDirections(2,TotalLegsBottom) = ChargeDirectionsUsedBottom(TotalLegsBottom);
                    end
                    
                    for jj = 2:(TotalLegsBottom-1)
                        NewStructure(1,TotalLegsBottom-jj+2) = TotalLegsBottom-jj;
                        NewChargeDirections(1,TotalLegsBottom-jj+2) = -1;
                        NewStructure(2,TotalLegsBottom-jj+1) = -BottomLegsOriginal(end-jj+1);
                        NewChargeDirections(2,TotalLegsBottom-jj+1) = ChargeDirectionsUsedBottom(TotalLegsBottom-jj+1);
                    end
                    
                    if TotalLegsBottom>=2
                        NewStructure(1,2) = -BottomLegsOriginal(1);
                        NewChargeDirections(1,2) = ChargeDirectionsUsedBottom(1);
                    end
                    
                    if TotalLegsTop == 1
                        NewStructure(2,1) = -TopLegsOriginal(1);
                        NewChargeDirections(2,1) = ChargeDirectionsUsedTop(1);
                    else
                        NewStructure(2,1) = TotalLegsBottom+TotalLegsTop-2;
                        NewStructure(2,TotalLegsBottom+TotalLegsTop-1) = -TopLegsOriginal(end);
                        NewChargeDirections(2,TotalLegsBottom+TotalLegsTop-1) = ChargeDirectionsUsedTop(TotalLegsTop);
                    end
                    
                    for jj = 2:(TotalLegsTop-1)
                        NewStructure(1,TotalLegsTop-jj+2+TotalLegsBottom-1) = TotalLegsTop-jj-1+TotalLegsBottom;
                        NewStructure(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = -TopLegsOriginal(end-jj+1);
                        NewChargeDirections(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = ChargeDirectionsUsedTop(TotalLegsTop-jj+1);
                    end
                    
                    if TotalLegsTop>=2
                        NewStructure(1,TotalLegsBottom+1) = -TopLegsOriginal(1);
                        NewChargeDirections(1,TotalLegsBottom+1) = ChargeDirectionsUsedTop(1);
                    end
                    
                   else
                    
                    if TotalLegsBottom == 0
                        
                        
                        
                    else %if TotalLegsTop == 0
                        
                    end
                
                end
               
               
            else
                NewStructure = zeros([2,0]);
                NewChargeDirections = zeros([2,0]);
            end
            
            if TotalLegsBottom>=2
               NewChargeSideInt = -ones([1,TotalLegsBottom-1]);
            else
                NewChargeSideInt = -ones([1,0]);
            end
            
            if TotalLegsTop>=2
               NewChargeSideInt = [NewChargeSideInt,ones([1,TotalLegsTop-1])];
            end
            
            if RenameLegs
                NewOrder = [BottomLegsOriginal,TopLegsOriginal];
                NewNewStructure = NewStructure;
                for kk = 1:numel(NewOrder)
                    NewNewStructure(NewStructure == -NewOrder(kk)) = -kk;
                end
                NewStructure = NewNewStructure;
            else
                NewOrder = 1:numel(self.CapType);
            end
            
            NewChargeSide = self.ChargeSide(NewOrder);
            NewCapType = self.CapType(NewOrder);
            
            
            self = self.SymReshape(NewStructure, NewChargeDirections, NewOrder, NewCapType, NewChargeSide, NewChargeSideInt, [], [], [], [], []);
        end
        
        function self = Braid(self, Braids,BraidDirections, RenameLegs)
            
            if nargin<4
                RenameLegs = true;
            end
            
            %note that this braids after all the other lined up braids and
            %replaces them, however we swap entries in CylicLegs rather
            %then CyclicLegs braided.
            
            AllNumber = numel(self.ChargeSide);
            
            NewStructure = self.Structure;
            NewChargeDirections = self.ChargeDirections;
            NewCapType = self.CapType;
            NewChargeSide = self.ChargeSide;
            NewChargeSideInt = self.ChargeSideInt;
            NewBraidings = self.Braidings;
            NewBraidingDirections = self.BraidingDirections;
            
            ExtraBraidings = Braids;
            ExtraBraidingDirections = BraidDirections;
            
            NewOrdering = 1:AllNumber;
            
            for kk = 1:numel(Braids)
                if Braids(kk) == AllNumber
                    ASite = self.CyclicLegsBraided(1);
                    BSite = self.CyclicLegsBraided(end);
                else
                    ASite = self.CyclicLegsBraided(Braids(kk)+1);
                    BSite = self.CyclicLegsBraided(Braids(kk));
                end
                New2Ordering = NewOrdering;
                New2Ordering([ASite,BSite]) = New2Ordering([BSite,ASite]);
                %New2Ordering(NewOrdering == BSite) = New2Ordering(NewOrdering == ASite);
                NewOrdering = New2Ordering;
                
                ASite = find(NewStructure(:)==-ASite);
                BSite = find(NewStructure(:)==-BSite);
                
                
                if numel(ASite)~=1
                    error('Affirmation Error: ASite should have exactly one entry');
                end
                if numel(BSite)~=1
                    error('Affirmation Error: BSite should have exactly one entry');
                end
                
                New2ChargeDirections = NewChargeDirections;
                New2ChargeDirections([ASite,BSite]) = New2ChargeDirections([BSite,ASite]);
                %New2Ordering(NewOrdering == BSite) = New2Ordering(NewOrdering == ASite);
                NewChargeDirections = New2ChargeDirections;
                
                %New2ChargeDirections = NewChargeDirections;
                %New2ChargeDirections(NewStructure == -ASite) = NewChargeDirections(NewStructure == -BSite);
                %New2ChargeDirections(NewStructure == -BSite) = NewChargeDirections(NewStructure == -ASite);
                %NewChargeDirections = New2ChargeDirections;
            end
            
            if ~RenameLegs
                NewNewStructure = NewStructure;
                for kk = 1:numel(NewOrdering)
                    NewNewStructure(NewStructure == -NewOrdering(kk)) = -kk;
                    %NewNewStructure(NewStructure == -NewOrdering(kk)) = New(NewStructure == -kk);
                end
                NewStructure = NewNewStructure;
                NewOrdering = 1:AllNumber;
            end
            
            self = self.SymReshape(NewStructure, NewChargeDirections, NewOrdering, NewCapType, NewChargeSide, NewChargeSideInt, NewBraidings, NewBraidingDirections, [], ExtraBraidings, ExtraBraidingDirections);
        end
        
        function self = RotateCW(self, Rotate, RenameLegs)
            
            %if Rotate is a cell then it should be the names of the indicies
            
            if nargin<3
                RenameLegs = false;
            end
            
            if iscell(Rotate)
                FlagRotateIsNames = true;
                Rotate = cell2mat(Rotate);
                %this means rotate is refering to the names rather the
                %locations
            else
                FlagRotateIsNames = false;
                %this means that rotate is refering to locations in the
                %cycle rather then the names
            end
            
            if ~isnumeric(Rotate)
                error('Error:Rotate should be numeric')
            end
            
            if ~all(IsInteger(Rotate(:))) || any(Rotate<=0) || any(Rotate>numel(self.ChargeSide))
                error('Error:Rotate should be positive integers less then or equal to the number of legs of the tensor')
            end
            
            
            Rotate = sort(Rotate(:),'ascend');
            if any(Rotate(1:(end-1))-Rotate(2:end)==0)
                error('Error: We can only rotate a leg once at a time');
            end
            
            
            Numbers = 1:numel(self.ChargeSide);
            TopLegsOriginal = Numbers(self.ChargeSide == +1)';
            TopLegsOriginal = self.CyclicLegs(any(repmat(self.CyclicLegs, [numel(TopLegsOriginal),1]) == repmat(TopLegsOriginal, [1,numel(self.CyclicLegs)]),1));
            
            BottomLegsOriginal = Numbers(self.ChargeSide == -1)';
            BottomLegsOriginal = self.CyclicLegs(any(repmat(self.CyclicLegs, [numel(BottomLegsOriginal),1]) == repmat(BottomLegsOriginal, [1,numel(self.CyclicLegs)]),1));
            
            
            if FlagRotateIsNames
                %This is if we signalled that the labels we will be using
                %are names rather then locations
                
                %now check that the rotate values given agree:
                RotateFromTop = any(repmat(TopLegsOriginal, [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(TopLegsOriginal)]),1);
                RotateFromBottom = any(repmat(BottomLegsOriginal, [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(BottomLegsOriginal)]),1);
                
            else
                %This is if we signalled that the labels we will be using
                %are locations rather then names
                
                %now check that the rotate values given agree:
                RotateFromTop = any(repmat(1:numel(TopLegsOriginal), [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(TopLegsOriginal)]),1);
                RotateFromBottom = any(repmat(numel(TopLegsOriginal)+(1:numel(BottomLegsOriginal)), [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(BottomLegsOriginal)]),1);
                
            end
            %this give the original values, now need to check that they
            %agree with my thoughts on this region (phases will be absorbed
            %by CCW)
            
            
            if any(RotateFromTop(1:(end-1))>RotateFromTop(2:end)) || any(RotateFromBottom(1:(end-1))>RotateFromBottom(2:end))
                error('Error: Rotate should only consider a series of legs that can be rotated clockwise');
            end
            
            AmountRotated = sum(RotateFromBottom);
            
            BottomLegsRotated = BottomLegsOriginal(RotateFromBottom);
            BottomLegsOriginal(RotateFromBottom) = [];
            
            TopLegsRotated = TopLegsOriginal(RotateFromTop);
            TopLegsOriginal(RotateFromTop) = [];
            
            TopLegsOriginal = [BottomLegsRotated,TopLegsOriginal];
            BottomLegsOriginal = [TopLegsRotated,BottomLegsOriginal];
            
            RotateNames = [TopLegsRotated, BottomLegsRotated];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %now we have the ordered structure and the rest is like
            %BringToRegular:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            CyclicLegsBraidedOriginal = [TopLegsOriginal, BottomLegsOriginal];
            %this is the order that I need the new legs to be in
            
            BottomLegsOriginal = BottomLegsOriginal(end:-1:1);
            
            ExternalChargeDirections = zeros([numel(CyclicLegsBraidedOriginal),1]);
            for kk = 1:numel(CyclicLegsBraidedOriginal)
                ExternalChargeDirections(kk) = self.ChargeDirections(self.Structure == -kk);
            end
            
            TotalLegsBottom = numel(BottomLegsOriginal);
            TotalLegsTop = numel(TopLegsOriginal);
            ChargeDirectionsUsedBottom = ExternalChargeDirections(BottomLegsOriginal);
            ChargeDirectionsUsedTop = ExternalChargeDirections(TopLegsOriginal);
            
            if numel(CyclicLegsBraidedOriginal)>=2
               NewStructure = zeros([2,numel(CyclicLegsBraidedOriginal)-1]);
               NewChargeDirections = ones([2,numel(CyclicLegsBraidedOriginal)-1]);
               if TotalLegsBottom*TotalLegsTop >0
                    NewStructure = zeros([2,TotalLegsBottom+TotalLegsTop-1]);
                    NewChargeDirections = ones([2,TotalLegsBottom+TotalLegsTop-1]);
                    
                    if TotalLegsBottom == 1
                        NewStructure(1,1) = -BottomLegsOriginal(1);
                        NewChargeDirections(1,1) = ChargeDirectionsUsedBottom(1);
                    else
                        NewStructure(1,1) = TotalLegsBottom-1;
                        NewChargeDirections(1,1) = -1;
                        NewStructure(2,TotalLegsBottom) = -BottomLegsOriginal(end);
                        NewChargeDirections(2,TotalLegsBottom) = ChargeDirectionsUsedBottom(TotalLegsBottom);
                    end
                    
                    for jj = 2:(TotalLegsBottom-1)
                        NewStructure(1,TotalLegsBottom-jj+2) = TotalLegsBottom-jj;
                        NewChargeDirections(1,TotalLegsBottom-jj+2) = -1;
                        NewStructure(2,TotalLegsBottom-jj+1) = -BottomLegsOriginal(end-jj+1);
                        NewChargeDirections(2,TotalLegsBottom-jj+1) = ChargeDirectionsUsedBottom(TotalLegsBottom-jj+1);
                    end
                    
                    if TotalLegsBottom>=2
                        NewStructure(1,2) = -BottomLegsOriginal(1);
                        NewChargeDirections(1,2) = ChargeDirectionsUsedBottom(1);
                    end
                    
                    if TotalLegsTop == 1
                        NewStructure(2,1) = -TopLegsOriginal(1);
                        NewChargeDirections(2,1) = ChargeDirectionsUsedTop(1);
                    else
                        NewStructure(2,1) = TotalLegsBottom+TotalLegsTop-2;
                        NewStructure(2,TotalLegsBottom+TotalLegsTop-1) = -TopLegsOriginal(end);
                        NewChargeDirections(2,TotalLegsBottom+TotalLegsTop-1) = ChargeDirectionsUsedTop(TotalLegsTop);
                    end
                    
                    for jj = 2:(TotalLegsTop-1)
                        NewStructure(1,TotalLegsTop-jj+2+TotalLegsBottom-1) = TotalLegsTop-jj-1+TotalLegsBottom;
                        NewStructure(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = -TopLegsOriginal(end-jj+1);
                        NewChargeDirections(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = ChargeDirectionsUsedTop(TotalLegsTop-jj+1);
                    end
                    
                    if TotalLegsTop>=2
                        NewStructure(1,TotalLegsBottom+1) = -TopLegsOriginal(1);
                        NewChargeDirections(1,TotalLegsBottom+1) = ChargeDirectionsUsedTop(1);
                    end
                    
                    NewChargeSideInt = ones([1,0]);
                    if TotalLegsBottom>=2
                       NewChargeSideInt = -ones([1,TotalLegsBottom-1]);
                    end
                    
                    if TotalLegsTop>=2
                       NewChargeSideInt = [NewChargeSideInt,ones([1,TotalLegsTop-1])];
                    end
                    
               else
                    
                    
                    if TotalLegsBottom == 0
                        NewChargeSideInt = ones([1,TotalLegsTop-2]);
                        
                        if TotalLegsTop == 2
                            NewStructure(1,1) = -TopLegsOriginal(1);
                            NewChargeDirections(1,1) = ChargeDirectionsUsedTop(1);
                            NewStructure(2,1) = -TopLegsOriginal(2);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedTop(2);
                        else
                            NewStructure(1,1) = TotalLegsTop-2;
                            NewStructure(2,1) = -TopLegsOriginal(end);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedTop(end);
                            NewStructure(1,2) = -TopLegsOriginal(1);
                            NewChargeDirections(1,2) = ChargeDirectionsUsedTop(1);
                            NewStructure(2,2) = -TopLegsOriginal(2);
                            NewChargeDirections(2,2) = ChargeDirectionsUsedTop(2);
                        end
                         
                        for jj = 2:(TotalLegsTop-2)
                            NewStructure(1,TotalLegsTop-jj+1) = TotalLegsTop-jj-1;
                            NewStructure(2,TotalLegsTop-jj+1) = -TopLegsOriginal(end-jj+1);
                            NewChargeDirections(2,TotalLegsTop-jj+1) = ChargeDirectionsUsedTop(TotalLegsTop-jj+1);
                        end
                        
                    else %if TotalLegsTop == 0
                        NewChargeSideInt = -ones([1,TotalLegsBottom-2]);
                        
                        if TotalLegsBottom == 2
                            NewStructure(1,1) = -BottomLegsOriginal(1);
                            NewChargeDirections(1,1) = ChargeDirectionsUsedBottom(1);
                            NewStructure(2,1) = -BottomLegsOriginal(2);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedBottom(2);
                        else
                            NewStructure(1,1) = TotalLegsBottom-2;
                            NewChargeDirections(1,1) = -1;
                            NewStructure(2,1) = -BottomLegsOriginal(end);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedBottom(end);
                            NewStructure(1,2) = -BottomLegsOriginal(1);
                            NewChargeDirections(1,2) = ChargeDirectionsUsedBottom(1);
                            NewStructure(2,2) = -BottomLegsOriginal(2);
                            NewChargeDirections(2,2) = ChargeDirectionsUsedBottom(2);
                        end
                         
                        for jj = 2:(TotalLegsBottom-2)
                            NewStructure(1,TotalLegsBottom-jj+1) = TotalLegsBottom-jj-1;
                            NewChargeDirections(1,TotalLegsBottom-jj+1) = -1;
                            NewStructure(2,TotalLegsBottom-jj+1) = -BottomLegsOriginal(end-jj+1);
                            NewChargeDirections(2,TotalLegsBottom-jj+1) = ChargeDirectionsUsedBottom(TotalLegsBottom-jj+1);
                        end
                        
                    end
                end
               
               
            else
                NewStructure = zeros([2,0]);
                NewChargeDirections = zeros([2,0]);
                NewChargeSideInt = zeros([1,0]);
            end
            
            %self.CyclicLegsBraided = [BottomLegsOriginal,TopLegsOriginal];
            %self.CyclicLegs = [BottomLegsOriginal,TopLegsOriginal];
            
            if RenameLegs
                NewOrder = [BottomLegsOriginal,TopLegsOriginal];
                NewNewStructure = NewStructure;
                for kk = 1:numel(NewOrder)
                    NewNewStructure(NewStructure == -NewOrder(kk)) = -kk;
                end
                NewStructure = NewNewStructure;
            else
                NewOrder = 1:numel(self.CapType);
            end
            NewChargeSide = self.ChargeSide;
            NewChargeSide(RotateNames) = -NewChargeSide(RotateNames);
            
            NewChargeSide = NewChargeSide(NewOrder);
            NewCapType = self.CapType(NewOrder);
            
            NewBraidings = 1+mod(self.Braidings+AmountRotated-1,numel(NewChargeSide));
            NewBraidingDirections = self.BraidingDirections;
            
            
            
            self = self.SymReshape(NewStructure, NewChargeDirections, NewOrder, NewCapType, NewChargeSide, NewChargeSideInt, NewBraidings, NewBraidingDirections, [], [], [],+1);
        end
        
        function self = RotateCCW(self,Rotate, RenameLegs)
           
            %if Rotate is a cell then it should be the names of the indicies
            
            if nargin<3
                RenameLegs = false;
            end
            
            if iscell(Rotate)
                FlagRotateIsNames = true;
                Rotate = cell2mat(Rotate);
                %this means rotate is refering to the names rather the
                %locations
            else
                FlagRotateIsNames = false;
                %this means that rotate is refering to locations in the
                %cycle rather then the names
            end
            
            if ~isnumeric(Rotate)
                error('Error:Rotate should be numeric')
            end
            
            if ~all(IsInteger(Rotate(:))) || any(Rotate<=0) || any(Rotate>numel(self.ChargeSide))
                error('Error:Rotate should be positive integers less then or equal to the number of legs of the tensor')
            end
            
            
            Rotate = sort(Rotate(:),'ascend');
            if any(Rotate(1:(end-1))-Rotate(2:end)==0)
                error('Error: We can only rotate a leg once at a time');
            end
            
            
            Numbers = 1:numel(self.ChargeSide);
            TopLegsOriginal = Numbers(self.ChargeSide == +1)';
            TopLegsOriginal = self.CyclicLegs(any(repmat(self.CyclicLegs, [numel(TopLegsOriginal),1]) == repmat(TopLegsOriginal, [1,numel(self.CyclicLegs)]),1));
            
            BottomLegsOriginal = Numbers(self.ChargeSide == -1)';
            BottomLegsOriginal = self.CyclicLegs(any(repmat(self.CyclicLegs, [numel(BottomLegsOriginal),1]) == repmat(BottomLegsOriginal, [1,numel(self.CyclicLegs)]),1));
            
            
            if FlagRotateIsNames
                %This is if we signalled that the labels we will be using
                %are names rather then locations
                
                %now check that the rotate values given agree:
                RotateFromTop = any(repmat(TopLegsOriginal, [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(TopLegsOriginal)]),1);
                RotateFromBottom = any(repmat(BottomLegsOriginal, [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(BottomLegsOriginal)]),1);
                
            else
                %This is if we signalled that the labels we will be using
                %are locations rather then names
                
                %now check that the rotate values given agree:
                RotateFromTop = any(repmat(1:numel(TopLegsOriginal), [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(TopLegsOriginal)]),1);
                RotateFromBottom = any(repmat(numel(TopLegsOriginal)+(1:numel(BottomLegsOriginal)), [numel(Rotate),1]) == repmat(Rotate(:), [1,numel(BottomLegsOriginal)]),1);
                
            end
            %this give the original values, now need to check that they
            %agree with my thoughts on this region (phases will be absorbed
            %by CCW)
            
            
            if any(RotateFromTop(1:(end-1))<RotateFromTop(2:end)) || any(RotateFromBottom(1:(end-1))<RotateFromBottom(2:end))
                error('Error: Rotate should only consider a series of legs that can be rotated Counterclockwise');
            end
            
            AmountRotated = sum(RotateFromTop);
            
            BottomLegsRotated = BottomLegsOriginal(RotateFromBottom);
            BottomLegsOriginal(RotateFromBottom) = [];
            
            TopLegsRotated = TopLegsOriginal(RotateFromTop);
            TopLegsOriginal(RotateFromTop) = [];
            
            TopLegsOriginal = [TopLegsOriginal,BottomLegsRotated];
            BottomLegsOriginal = [BottomLegsOriginal,TopLegsRotated];
            
            RotateNames = [TopLegsRotated, BottomLegsRotated];
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %now we have the ordered structure and the rest is like
            %BringToRegular:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            CyclicLegsBraidedOriginal = [TopLegsOriginal, BottomLegsOriginal];
            %this is the order that I need the new legs to be in
            
            BottomLegsOriginal = BottomLegsOriginal(end:-1:1);
            
            ExternalChargeDirections = zeros([numel(CyclicLegsBraidedOriginal),1]);
            for kk = 1:numel(CyclicLegsBraidedOriginal)
                ExternalChargeDirections(kk) = self.ChargeDirections(self.Structure == -kk);
            end
            
            TotalLegsBottom = numel(BottomLegsOriginal);
            TotalLegsTop = numel(TopLegsOriginal);
            ChargeDirectionsUsedBottom = ExternalChargeDirections(BottomLegsOriginal);
            ChargeDirectionsUsedTop = ExternalChargeDirections(TopLegsOriginal);
            
            if numel(CyclicLegsBraidedOriginal)>=2
               NewStructure = zeros([2,numel(CyclicLegsBraidedOriginal)-1]);
               NewChargeDirections = ones([2,numel(CyclicLegsBraidedOriginal)-1]);
               
               if TotalLegsBottom*TotalLegsTop >0
                    
                    NewChargeSideInt = [-ones([1,TotalLegsBottom-1]),ones([1,TotalLegsTop-1])];
                    if TotalLegsBottom == 1
                        NewStructure(1,1) = -BottomLegsOriginal(1);
                        NewChargeDirections(1,1) = ChargeDirectionsUsedBottom(1);
                    else
                        NewStructure(1,1) = TotalLegsBottom-1;
                        NewChargeDirections(1,1) = -1;
                        NewStructure(2,TotalLegsBottom) = -BottomLegsOriginal(end);
                        NewChargeDirections(2,TotalLegsBottom) = ChargeDirectionsUsedBottom(end);
                    end
                    
                    for jj = 2:(TotalLegsBottom-1)
                        NewStructure(1,TotalLegsBottom-jj+2) = TotalLegsBottom-jj;
                        NewChargeDirections(1,TotalLegsBottom-jj+2) = -1;
                        NewStructure(2,TotalLegsBottom-jj+1) = -BottomLegsOriginal(end-jj+1);
                        NewChargeDirections(2,TotalLegsBottom-jj+1) = ChargeDirectionsUsedBottom(TotalLegsBottom-jj+1);
                    end
                    
                    if TotalLegsBottom>=2
                        NewStructure(1,2) = -BottomLegsOriginal(1);
                        NewChargeDirections(1,2) = ChargeDirectionsUsedBottom(1);
                    end
                    
                    if TotalLegsTop == 1
                        NewStructure(2,1) = -TopLegsOriginal(1);
                        NewChargeDirections(2,1) = ChargeDirectionsUsedTop(1);
                    else
                        NewStructure(2,1) = TotalLegsBottom+TotalLegsTop-2;
                        NewStructure(2,TotalLegsBottom+TotalLegsTop-1) = -TopLegsOriginal(end);
                        NewChargeDirections(2,TotalLegsBottom+TotalLegsTop-1) = ChargeDirectionsUsedTop(end);
                    end
                    
                    for jj = 2:(TotalLegsTop-1)
                        NewStructure(1,TotalLegsTop-jj+2+TotalLegsBottom-1) = TotalLegsTop-jj-1+TotalLegsBottom;
                        NewStructure(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = -TopLegsOriginal(end-jj+1);
                        NewChargeDirections(2,TotalLegsTop-jj+1+TotalLegsBottom-1) = ChargeDirectionsUsedTop(TotalLegsTop-jj+1);
                    end
                    
                    if TotalLegsTop>=2
                        NewStructure(1,TotalLegsBottom+1) = -TopLegsOriginal(1);
                        NewChargeDirections(1,TotalLegsBottom+1) = ChargeDirectionsUsedTop(1);
                    end
                    
                    
                    NewChargeSideInt = ones([1,0]);
                    if TotalLegsBottom>=2
                       NewChargeSideInt = -ones([1,TotalLegsBottom-1]);
                    end
                    
                    if TotalLegsTop>=2
                       NewChargeSideInt = [NewChargeSideInt,ones([1,TotalLegsTop-1])];
                    end
                    
               else
                    
                    if TotalLegsBottom == 0
                        NewChargeSideInt = ones([1,TotalLegsTop-2]);
                        
                        if TotalLegsTop == 2
                            NewStructure(1,1) = -TopLegsOriginal(1);
                            NewChargeDirections(1,1) = ChargeDirectionsUsedTop(1);
                            NewStructure(2,1) = -TopLegsOriginal(2);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedTop(2);
                        else
                            NewStructure(1,1) = TotalLegsTop-2;
                            NewStructure(2,1) = -TopLegsOriginal(end);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedTop(end);
                            NewStructure(1,2) = -TopLegsOriginal(1);
                            NewChargeDirections(1,2) = ChargeDirectionsUsedTop(1);
                            NewStructure(2,2) = -TopLegsOriginal(2);
                            NewChargeDirections(2,2) = ChargeDirectionsUsedTop(2);
                        end
                         
                        for jj = 2:(TotalLegsTop-2)
                            NewStructure(1,TotalLegsTop-jj+1) = TotalLegsTop-jj-1;
                            NewStructure(2,TotalLegsTop-jj+1) = -TopLegsOriginal(end-jj+1);
                            NewChargeDirections(2,TotalLegsTop-jj+1) = ChargeDirectionsUsedTop(TotalLegsTop-jj+1);
                        end
                        
                    else %if TotalLegsTop == 0
                        NewChargeSideInt = -ones([1,TotalLegsBottom-2]);
                        
                        if TotalLegsBottom == 2
                            NewStructure(1,1) = -BottomLegsOriginal(1);
                            NewChargeDirections(1,1) = ChargeDirectionsUsedBottom(1);
                            NewStructure(2,1) = -BottomLegsOriginal(2);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedBottom(2);
                        else
                            NewStructure(1,1) = TotalLegsBottom-2;
                            NewChargeDirections(1,1) = -1;
                            NewStructure(2,1) = -BottomLegsOriginal(end);
                            NewChargeDirections(2,1) = ChargeDirectionsUsedBottom(end);
                            NewStructure(1,2) = -BottomLegsOriginal(1);
                            NewChargeDirections(1,2) = ChargeDirectionsUsedBottom(1);
                            NewStructure(2,2) = -BottomLegsOriginal(2);
                            NewChargeDirections(2,2) = ChargeDirectionsUsedBottom(2);
                        end
                         
                        for jj = 2:(TotalLegsBottom-2)
                            NewStructure(1,TotalLegsBottom-jj+1) = TotalLegsBottom-jj-1;
                            NewChargeDirections(1,TotalLegsBottom-jj+1) = -1;
                            NewStructure(2,TotalLegsBottom-jj+1) = -BottomLegsOriginal(end-jj+1);
                            NewChargeDirections(2,TotalLegsBottom-jj+1) = ChargeDirectionsUsedBottom(TotalLegsBottom-jj+1);
                        end
                        
                    end
                
               end
                
               
            else
                NewStructure = zeros([2,0]);
                NewChargeDirections = zeros([2,0]);
                NewChargeSideInt = -ones([1,0]);
            end
            
            
            %self.CyclicLegsBraided = [BottomLegsOriginal,TopLegsOriginal];
            %self.CyclicLegs = [BottomLegsOriginal,TopLegsOriginal];
            
            if RenameLegs
                NewOrder = [BottomLegsOriginal,TopLegsOriginal];
                NewNewStructure = NewStructure;
                for kk = 1:numel(NewOrder)
                    NewNewStructure(NewStructure == -NewOrder(kk)) = -kk;
                end
                NewStructure = NewNewStructure;
            else
                NewOrder = 1:numel(self.CapType);
            end
            NewChargeSide = self.ChargeSide;
            NewChargeSide(RotateNames) = -NewChargeSide(RotateNames);
            
            NewChargeSide = NewChargeSide(NewOrder);
            NewCapType = self.CapType(NewOrder);
            
            NewBraidings = 1+mod(self.Braidings-AmountRotated-1,numel(NewChargeSide));
            NewBraidingDirections = self.BraidingDirections;
            
            
            
            self = self.SymReshape(NewStructure, NewChargeDirections, NewOrder, NewCapType, NewChargeSide, NewChargeSideInt, NewBraidings, NewBraidingDirections, [], [], [],-1);%,-AmountRotated);
        end
        
        function self = ModifyPhase(self,Legs, Directions)
            
            %this is the name of the leg not its location on the loop
            
            ModifyPhases = [Legs(:),Directions(:)];
            
            self = self.SymReshape(self.Structure, self.ChargeDirections, 1:numel(self.ChargeSide), self.CapType, self.ChargeSide, self.ChargeSideInt, ...
                self.Braidings, self.BraidingDirections, ModifyPhases, [], []);
            
        end
        
        function self = CombineLegs(self,Legs)
            
        end
        
        function self = EvaluateMoves(self, Moves)
            %this function evalues a set of moves on this tensor
            
            %if nargin<2; Moves = []; end
            
            if ~isempty(Moves)
                
                LabelsIn = Moves{1};
                LabelsOut = Moves{2};
                
                LabelsWorking = [self.ChargeLabelsExternal,self.ChargeLabelsInternal,self.MultiplicitiesInternal];
                
                [~,~,NumberingIn] = unique([LabelsIn;LabelsWorking],'rows','stable');
                
                if max(NumberingIn)>size(LabelsIn,1)
                    error('Error: I don''t have all possible evaluations')
                end
                
                NumberingIn = NumberingIn((size(LabelsIn,1)+1):end);
                
                NumbersInUsing = unique(NumberingIn);
                MovesUsing = Moves{3}(NumbersInUsing,:);
                
                [~,NumbersOutUsing,~] = find(MovesUsing);
                NumbersOutUsing = unique(NumbersOutUsing);
                
                LabelsOutUsing = LabelsOut(NumbersOutUsing,:);
                
                MovesUsing = MovesUsing(:,NumbersOutUsing);
                
                TensorsTemp = cell([numel(NumbersOutUsing),1]);
                TensorsIn = self.TensorEntries;
                
                for kk = 1:size(MovesUsing,1)
                    [~, Outs, OutsMult] = find(MovesUsing(kk,:));
                    WorkingTensor = TensorsIn{kk};
                    for aa = 1:numel(Outs)
                        jj = Outs;
                        
                        if isempty(TensorsTemp{jj})
                            TensorsTemp{jj} = OutsMult(aa)*WorkingTensor;
                        else
                            TensorsTemp{jj} = TensorsTemp{jj} + OutsMult(aa)*WorkingTensor;
                        end
                    end
                end
                
                if size(LabelsOutUsing,2)>2
                    N = ((size(LabelsOutUsing,2)+4)/3);
                    if ~IsInteger(N)
                        error('Affirmation Error: Something is wrong with the number of output Charges and Mults');
                    end
                    
                    self.ChargeLabelsExternal = LabelsOutUsing(:,1:N);
                    self.ChargeLabelsInternal = LabelsOutUsing(:,N+(1:(N-2)));
                    self.MultiplicitiesInternal = LabelsOutUsing(:,2*N-2+(1:(N-2)));
                else
                    self.ChargeLabelsExternal = LabelsOutUsing;
                    self.ChargeLabelsInternal = zeros([size(LabelsOutUsing,1),0]);
                    self.MultiplicitiesInternal = zeros([size(LabelsOutUsing,1),0]);
                end
                
                self.TensorEntries = TensorsTemp;
%                self.Structure = Moves{4};
%                self.ChargeDirections = Moves{5};
%                self.ChargeSide = Moves{6};
%                self.ChargeLabelsExternal
                
            end 
            
            %otherwise we do nothing to self and return it.
        end
        
        function self = FreePermute(self, Permute)
            %this function performs permutations on the tensor blocks
            
            %if nargin<2; Permute = []; end
                
            if ~isempty(Permute)
                
                if ~isnumeric(Permute)
                    error('Error: Input of FreePermute must be a number')
                end
                Permute = Permute(:);
                
                if any(~IsInteger(Permute))
                    error('Error: Input of FreePermute must be an integer')
                end
                
                if isequal(sort(Permute,'ascend'), (1:size(self.ChargeLabelsExternal,2))')
                    
                    self.TensorEntriesSizes = self.TensorEntriesSizes(:,Permute);
                    %self.FreePermuteTracker = self.FreePermuteTracker(Permute);
                    self.ChargeLegDimensions = self.ChargeLegDimensions(:,Permute);
                    
                    if numel(Permute) == 1
                        Permute = [1,2];
                    end
                    
                    for kk = 1:numel(self.TensorEntries)
                        self.TensorEntries{kk} = permute(self.TensorEntries{kk},Permute);
                    end
                    
                else
                    error('Error: Input of FreePermute must be a list of numbers from 1 to the number of outgoing legs of this tensor. Which in this case is %d', size(self.ChargeLabelsExternal,2))
                end
                
            end
            %otherwise we do nothing to self and return it.
        end
        
        function self = MakeAShape(self)
            
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{kk} = [];
            end
            self.FlagIsAShape = true;
            
        end
        
        function Sym = Symmetry(self)
            Sym = self.SymHandle;
        end
        
        function [Tensors, Sizes] = getTensorEntries(self)
            Tensors = self.TensorEntries;
            Sizes = self.TensorEntriesSizes;
        end
        
        function Legs = nLegs(self)
            Legs = size(self.ChargeLegDimensions,2);
        end
        
        function BoolVal = ConnectableLegs(Tensor1, Tensor2, Leg1, Leg2)
            BoolVal = true;
            Charges1 = Tensor1.ChargeAccess(Leg1);
            Charges2 = Tensor2.ChargeAccess(Leg2);
            if ~isequal(size(Charges1),size(Charges2))
                BoolVal = false;
                return;
            elseif ~all(Charges1 == Charges2 | Charges1 == 0 | Charges2 == 0)
                BoolVal = false;
                return;
            end
            
            if ~isequal(Tensor1.SymHandle,Tensor2.SymHandle)
                BoolVal = false;
                return;
            end
            
            if ~Tensor1.SymHandle.IsSelfDual
                %then direction does matter
                if Tensor1.ChargeDirections(Tensor1.Structure == -Leg1)~=-Tensor2.ChargeDirections(Tensor2.Structure == -Leg2)
                    BoolVal = false;
                end
            end
        end
        
        function Charges = ChargeAccess(self,LegNumber)
            if nargin<2
                Charges = self.ChargeLegDimensions;
            else
                if IsInteger(LegNumber)
                    LegNumber = round(LegNumber);
                    if (LegNumber >size(self.ChargeLegDimensions,2)) || LegNumber<=0
                        error('LegNumber must be between 1 and the number of legs on the tensor')
                    end
                    
                else
                    error('LegNumber must be an integer')
                end
                Charges = self.ChargeLegDimensions(:,LegNumber);
            end
        end
        
        function Out = minus(A,B)
            Out = A+((-1)*B);
        end
        
        function Out = plus(A,B)
            IsTensorA = isa(A,'SymTensor');
            IsTensorB = isa(B,'SymTensor');
            
            if IsTensorA
                if IsTensorB
                    %now both are tensors, check they have the same shape:
                    if (~isequal(size(A.ChargeLegDimensions),size(B.ChargeLegDimensions))||~isequal(A.SymHandle,B.SymHandle))
                        error('error: the tensors must be the same shape if you plan to add them')
                    elseif ~all(all(A.ChargeLegDimensions==B.ChargeLegDimensions|A.ChargeLegDimensions==0|B.ChargeLegDimensions==0))
                        error('error: the tensors must be the same shape if you plan to add them')
                    end
                    
                    A.ChargeLegDimensions(A.ChargeLegDimensions == 0) = B.ChargeLegDimensions(A.ChargeLegDimensions == 0);
                    
                    %now need to modify so that they have the same
                    %structure, ChargeLabels and ChargeDirections
                    %at the moment it is assumed that they must be the same
                    if (~isequal(A.Structure,B.Structure)||~isequal(A.ChargeDirections,B.ChargeDirections))
                        error('error: the tensors must be the same shape if you plan to add them')
                    end
                    
                    ChargesA = [A.ChargeLabelsExternal,A.ChargeLabelsInternal,A.MultiplicitiesInternal];
                    ChargesB = [B.ChargeLabelsExternal,B.ChargeLabelsInternal,B.MultiplicitiesInternal];
                    
                    [OutCharges,SizesIndex,Lists] = unique([ChargesA;ChargesB],'rows');
                    ListsA = Lists(1:size(ChargesA,1),:);
                    ListsB = Lists((size(ChargesA,1)+1):end,:);
                    
                    TensorEntriesSizesIn = [A.TensorEntriesSizes;B.TensorEntriesSizes];
                    
                    TensorEntriesSizesOut = TensorEntriesSizesIn(SizesIndex,:);
                    
                    OutChargeLabelsExternal = OutCharges(:,1:size(A.ChargeLabelsExternal,2));
                    OutChargeLabelsInternal = OutCharges(:,size(A.ChargeLabelsExternal,2)+(1:size(A.ChargeLabelsInternal,2)));
                    OutMultiplicitiesInternal = OutCharges(:,size(A.ChargeLabelsExternal,2)+size(A.ChargeLabelsInternal,2)+(1:size(A.MultiplicitiesInternal,2)));
                    
                    TensorEntriesOut = cell([size(OutChargeLabelsExternal,1),1]);
                    for kk = 1:size(ListsA,1)
                        TensorEntriesOut{ListsA(kk)} = A.TensorEntries{kk};
                    end
                    for kk = 1:size(ListsB,1)
                        if ~isempty(TensorEntriesOut{ListsB(kk)})
                            TensorEntriesOut{ListsB(kk)} = TensorEntriesOut{ListsB(kk)}+B.TensorEntries{kk};
                        else
                            TensorEntriesOut{ListsB(kk)} = B.TensorEntries{kk};
                        end
                    end
                    
                    A.ChargeLabelsExternal = OutChargeLabelsExternal;
                    A.ChargeLabelsInternal = OutChargeLabelsInternal;
                    A.MultiplicitiesInternal = OutMultiplicitiesInternal;
                    
                    A.TensorEntries = TensorEntriesOut;
                    A.TensorEntriesSizes = TensorEntriesSizesOut;
                    Out = A;
                else
                    if ~isnumeric(B)
                        error('error: B must be a number')
                    end
                    if numel(B)~=1
                        error('error: A must be a single number')
                    end
                    for kk = 1:numel(A.TensorEntries)
                        A.TensorEntries{kk} = A.TensorEntries{kk}+B;
                    end
                    Out = A;
                end
            else
                if ~isnumeric(A)
                    error('error: A must be a number')
                end
                if numel(A)~=1
                    error('error: A must be a single number')
                end
                for kk = 1:numel(B.TensorEntries)
                    B.TensorEntries{kk} = B.TensorEntries{kk}+A;
                end
                Out = B;
            end
        end
        
        function Out = times(A,B)
            IsTensorA = isa(A,'SymTensor');
            IsTensorB = isa(B,'SymTensor');
            
            if IsTensorA
                if IsTensorB
                    error('error: can''t multiply two tensors, you need to use SymCon for that')
                else
                    if ~isnumeric(B)
                        error('error: B must be a number')
                    end
                    if numel(B)~=1
                        error('error: A must be a single number')
                    end
                    for kk = 1:numel(A.TensorEntries)
                        A.TensorEntries{kk} = A.TensorEntries{kk}*B;
                    end
                    Out = A;
                end
            else
                if ~isnumeric(A)
                    error('error: A must be a number')
                end
                if numel(A)~=1
                    error('error: A must be a single number')
                end
                for kk = 1:numel(B.TensorEntries)
                    B.TensorEntries{kk} = B.TensorEntries{kk}*A;
                end
                Out = B;
            end
        end

        function Out = mtimes(A,B)
            IsTensorA = isa(A,'SymTensor');
            IsTensorB = isa(B,'SymTensor');
            
            if IsTensorA
                if IsTensorB
                    error('error: can''t multiply two tensors, you need to use SymCon for that')
                else
                    if ~isnumeric(B)
                        error('error: B must be a number')
                    end
                    if numel(B)~=1
                        error('error: B must be a single number')
                    end
                    for kk = 1:numel(A.TensorEntries)
                        A.TensorEntries{kk} = A.TensorEntries{kk}*B;
                    end
                    Out = A;
                end
            else
                if ~isnumeric(A)
                    error('error: A must be a number')
                end
                if numel(A)~=1
                    error('error: A must be a single number')
                end
                for kk = 1:numel(B.TensorEntries)
                    B.TensorEntries{kk} = B.TensorEntries{kk}*A;
                end
                Out = B;
            end
        end
        
        function Out = transpose(A)
            Out = A.TenConj(true, false);
        end
        
        
        function Out = ctranspose(A)
            Out = A.TenConj(true);
        end
        
        function Out = getNumber(self)
            if self.FlagIsANumber
                Out = self.TensorEntries{1}(1);
            else
                Out = NaN;
            end
        end
        
        function Out = getTotalDimension(self, LegNumbers)
            Out = prod(sum(self.ChargeLegDimensions(:,LegNumbers).*repmat(self.SymHandle.Dim(:,2),[1,numel(LegNumbers)]),1),2);
        end
        
        function Out = Trace(self, Legs, Shape)
            %if ~self.SymHandle.IsSelfDual
                InverseIrrep = self.SymHandle.getInverseIrrep;
                if isequal(self.ChargeSide, [-ones([1,Legs(1)]),ones([1,Legs(2)])])
                    FlagCorrectSides = true;
                elseif isequal(self.ChargeSide, [ones([1,Legs(1)]),-ones([1,Legs(2)])])
                    FlagCorrectSides = false;
                    self.ChargeSide = -self.ChargeSide;
                    self.ChargeLabelsInternal = InverseIrrep(self.ChargeLabelsInternal);
                    self.ChargeDirections(self.Structure>0) = -self.ChargeDirections(self.Structure>0);
                else
                    error('Error: The Tensor Trace doesn''t work unless all the legs of a side are grouped up at any one time')
                end
            %end
            if nargin<3
                [Matrix,Data] = SymTensor.SymTen2Mat(self, Legs);
            else
                [Matrix,Data] = SymTensor.SymTen2Mat(self, Legs,Shape);
            end
            
            if ~isempty(Data{8})
                error('Error: Tensor Trace doesn''t work with currently braided tensors, the unresolved braiding has to be removed before we can proceed');
            end
            
            Data = Data{1};
            
            Out = 0;
            for dd = 1:size(self.SymHandle.Dim,1)
                Loc = find(Data == dd);
                if ~isempty(Loc)
                    if numel(Loc)>1
                        error('Affirmation Error: Too many entries in Loc to use');
                    end
                    Out = Out+self.SymHandle.Dim(dd,2)*trace(Matrix{Loc});
                end
            end
            
        end
        
        function Out = SymCut(self, Accuracy)
            if nargin<2
                Accuracy = 10^-14;
            end
            
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{1} = (abs(self.TensorEntries{1})>Accuracy)*1;
            end
            Out = self;
            
        end
        
        function Out = SymMax(self)
            Out = 0;
            for kk = 1:numel(self.TensorEntries)
                Out = max(Out, max(self.TensorEntries{kk}(:)));
            end
        end
    end
    
    methods(Access = 'private')
        
        function self = SymTensor(TensorEntries, Structure, ChargeDirections, MultiplicitiesInternal, ChargeLabelsInternal,...
                ChargeLabelsExternal, Braidings, BraidingDirections,StoredLocations,ChargeLegDimensions,SymHandle,ChargeSide,ChargeSideInt, CapType)
            if nargin<14
                CapType = zeros([1,size(ChargeLabelsExternal,2)]);
            end
            self.TensorEntries = TensorEntries;
            self.Structure = Structure;
            self.ChargeLabelsInternal = ChargeLabelsInternal;
            self.ChargeLabelsExternal = ChargeLabelsExternal;
            self.MultiplicitiesInternal = MultiplicitiesInternal;
            self.Braidings = Braidings;
            self.BraidingDirections = BraidingDirections;
            self.ChargeDirections = ChargeDirections;
            self.ChargeLegDimensions = ChargeLegDimensions;
            self.SymHandle = SymHandle;
            self.ChargeSide = ChargeSide;
            self.ChargeSideInt = ChargeSideInt;
            
            self.CapType = CapType;
            
            Numbers = 1:numel(ChargeSide);
            self.CyclicLegsBraided = Numbers(ChargeSide == -1);
            self.CyclicLegsBraided = [Numbers(ChargeSide == +1), self.CyclicLegsBraided(end:-1:1)];
            
            TempNumbers = self.CyclicLegsBraided;
            if ~isempty(Braidings)
                for kk = 1:numel(Braidings)
                    ll = Braidings(kk);
                    if ll~=Numbers(end);
                        TempNumbers([ll+1,ll]) = TempNumbers([ll,ll+1]);
                    else
                        TempNumbers([1,end]) = TempNumbers([end,1]);
                    end
                end
            end
            self.CyclicLegs = TempNumbers;%self.CyclicLegsBraided(NumbersReverseBraiding);
            
            
            %if isempty(StoredLocations)
            %we always compute StoredLocations as then we can check that
            %the structure is allowed.
            for jj = max(Structure(:)):-1:1
                Updating = jj;
                Using = [];
                while ~isempty(Updating)
                    kk = Updating(end);
                    Updating(end) = [];
                    Use1 = Structure(1,kk+1);
                    Use2 = Structure(2,kk+1);
                    if Use1>0
                        if Use1>jj
                            Using = [Using,StoredLocations{Use1}];
                        else
                            Updating = [Updating,Use1];
                        end
                    else %Use1<0
                        Using = [Using, -Use1];
                    end
                    
                    if Use2>0
                        if Use2>jj
                            Using = [Using,StoredLocations{Use2}];
                        else
                            Updating = [Updating,Use2];
                        end
                    else %Use2<0
                        Using = [Using, -Use2];
                    end
                end
                StoredLocations{jj} = Using;
                
                
                %now checking, that these are all next to each other
                %note that Using has at least 2 elements
                
                WorkOn = Using(1);
                Using(1) = [];
                while ~isempty(Using)
                    CurrentUse = WorkOn(end);
                    
                    Temp = find(WorkOn(end) == self.CyclicLegsBraided,1,'first');
                    WorkOn(end) = [];
                    
                    if numel(Temp)~=1
                        error('Affirmation Error: We should always be able to find this')
                    end
                    Temp = mod(Temp+[0,-2],Numbers(end))+1;
                    
                    TempFind1 = find(Using==self.CyclicLegsBraided(Temp(1)));
                    
                    if numel(TempFind1)>1
                        error('Affirmation Error: We shouldn''t be able to find more then 1')
                    end
                    
                    if ~isempty(TempFind1)
                        WorkOn = [WorkOn,Using(TempFind1)];
                        Using(TempFind1) = [];
                    end
                    
                    TempFind2 = find(Using==self.CyclicLegsBraided(Temp(2)));
                    
                    if numel(TempFind2)>1
                        error('Affirmation Error: We shouldn''t be able to find more then 1')
                    end
                    
                    if ~isempty(TempFind2)
                        WorkOn = [WorkOn,Using(TempFind2)];
                        Using(TempFind2) = [];
                    end
                    
                    if isempty(WorkOn)
                        error('Error: This is not a permissable tensor');
                    end
                
                    
                end
                
                
                
            end
            
            self.StoredLocations = StoredLocations;
            
            if ~isempty(TensorEntries)
                TotalLegs = size(ChargeLabelsExternal,2);
                TensorEntriesSizes = ones([size(TensorEntries,1), TotalLegs]);
                for kk = 1:size(TensorEntries,1)
                    TempSizing = size(TensorEntries{kk});
                    if TotalLegs == 1;
                        TensorEntriesSizes(kk,1) = prod(TempSizing);
                    else
                        TensorEntriesSizes(kk,1:numel(TempSizing)) = TempSizing;
                    end
                end
                
                if TotalLegs == 0
                    TensorEntriesSizes = zeros([1,0]);
                end
            else
                TotalLegs = size(ChargeLabelsExternal,2);
                TensorEntriesSizes = ones([0, TotalLegs]);
                
                if TotalLegs == 0
                    TensorEntriesSizes = zeros([1,0]);
                end
            end
            
            
            
            self.TensorEntriesSizes = TensorEntriesSizes;
            
            
            FlagIsANumber = false;
            FlagIsNotSymmetric = true;
            if any(any(ChargeLabelsExternal~=SymHandle.getTrivialIrrep)) 
                %if all external charges are trivial then all internal ones are as well, so we don't need this part 
                %|| any(ChargeLabelsInternal~=SymHandle.getTrivialIrrep)
                FlagIsNotSymmetric = false;
            else
                if any(MultiplicitiesInternal ~=1)
                    error('There shouldn''t be more then one Multiplicity given that everything is trivial')
                end
                
                if numel(TensorEntries)>1
                    error('There shouldn''t be more then one tensor entry given that everthing is trivial')
                end
                
%                if numel(TensorEntries)==0
%                    error('This tensor is empty')
%                end
                
                %now check if only a single number
                if isempty(ChargeLegDimensions)
                    FlagIsANumber = true;
                elseif prod(ChargeLegDimensions) == 1
                    FlagIsANumber = true;
                end
                
                
            end
            
            self.FlagIsNotSymmetric = FlagIsNotSymmetric;
            
            self.FlagIsANumber = FlagIsANumber;
            
            %self.FreePermuteTracker = 1:TotalLegs;
            
        end
    end
    
    methods(Access = 'public', Static = true)
        
        function [ChargeLabelsExternalUsed, ChargeLabelsInternalUsed,MultiplicitiesInternalUsed] = CreateUniqueEntries(AllowedSingleCharges, ChargeDirectionsUsed, TotalLegs,FlagNonAbelian, FlagMultiplicities, Fusion2)
            
            
            %now work out the combinations allowed
            
            ChargeLabelsExternalUsed = AllowedSingleCharges(:,1);
            
            if TotalLegs>1
            
            for jj = 2:TotalLegs
                
                ChargeLabelsExternalUsed = [repmat(ChargeLabelsExternalUsed, [size(AllowedSingleCharges,1), 1]),...
                    reshape(repmat(AllowedSingleCharges(:,jj)', [size(ChargeLabelsExternalUsed,1),1]),[size(ChargeLabelsExternalUsed,1)*size(AllowedSingleCharges,1),1])];
                
            end
            
            ChargeLabelsInternalUsed = zeros([size(ChargeLabelsExternalUsed,1),size(ChargeLabelsExternalUsed,2)-1]);
            MultiplicitiesInternalUsed = ones([size(ChargeLabelsExternalUsed,1),size(ChargeLabelsExternalUsed,2)-1]);
            
            if ~FlagNonAbelian
                
                %in the abelian case we know that the sum is fixed as one
                %output, we want to store that as we are using the Symmetry
                %to fix the dimension, may add an additional option later
                %where we store if we have an abelian model even if it
                %isn't know, however this will be important in cases where
                %we have braiding where knowing the internal structure is a
                %must, so it would have to be mutually exclusive to the
                %symmetry containing braiding in general.
                                
                %first one is special:
                %FIXHERE need to take case of one leg each side into
                %account
                for jj = size(ChargeLabelsExternalUsed,1):-1:1
                    if ~any(ChargeDirectionsUsed([1,2]) ==-1)
                        
                        Temp = find(Fusion2(:,ChargeLabelsExternalUsed(jj,1),ChargeLabelsExternalUsed(jj,2))~=0);
                        
                    elseif ChargeDirectionsUsed(1)*ChargeDirectionsUsed(2) == -1
                        
                        if ChargeDirectionsUsed(1) == -1
                            Temp = find(Fusion2(:,InverseIrrep(ChargeLabelsExternalUsed(jj,1)),ChargeLabelsExternalUsed(jj,2))~=0);
                        else
                            Temp = find(Fusion2(:,ChargeLabelsExternalUsed(jj,1),InverseIrrep(ChargeLabelsExternalUsed(jj,2)))~=0);
                        end
                        
                    else
                        Temp = find(Fusion2(:,InverseIrrep(ChargeLabelsExternalUsed(jj,1)),InverseIrrep(ChargeLabelsExternalUsed(jj,2)))~=0);
                    end
                    if ~isempty(Temp)
                        ChargeLabelsInternalUsed(jj,1) = Temp;
                    else
                        ChargeLabelsInternalUsed(jj,:) = [];
                        ChargeLabelsExternalUsed(jj,:) = [];
                    end
                    
                end
                
                for kk = 1:(size(ChargeLabelsExternalUsed,2)-2)
                    
                    for jj = size(ChargeLabelsExternalUsed,1):-1:1
                        if ChargeDirectionsUsed(kk+2) ==1
                            Temp = find(Fusion2(:,ChargeLabelsInternalUsed(jj,kk),ChargeLabelsExternalUsed(jj,kk+2))~=0);
                        else
                            Temp = find(Fusion2(:,ChargeLabelsInternalUsed(jj,kk),InverseIrrep(ChargeLabelsExternalUsed(jj,kk+2)))~=0);
                        end
                        
                        if ~isempty(Temp)
                            ChargeLabelsInternalUsed(jj,kk+1) = Temp;
                        else
                            ChargeLabelsInternalUsed(jj,:) = [];
                            ChargeLabelsExternalUsed(jj,:) = [];
                        end
                        
                    end
                    
                end
                
                %Now we fix up the remainder
                
                MultiplicitiesInternalUsed = ones([size(ChargeLabelsInternalUsed,1), TotalLegs-2]);
                
                
            else
                
                %if we have a non-abelian model:
                
                %if we don't need to worry about multiplicities, then we
                %only have one copy of each Temp
                if ~FlagMultiplicities
                    %first one is special:
                    for jj = size(ChargeLabelsExternalUsed,1):-1:1
                        if ~any(ChargeDirectionsUsed([1,2]) ==-1)
                            
                            Temp = find(Fusion2(:,ChargeLabelsExternalUsed(jj,1),ChargeLabelsExternalUsed(jj,2))~=0);
                            
                        elseif ChargeDirectionsUsed(1)*ChargeDirectionsUsed(2) == -1
                            
                            if ChargeDirectionsUsed(1) == -1
                                Temp = find(Fusion2(:,InverseIrrep(ChargeLabelsExternalUsed(jj,1)),ChargeLabelsExternalUsed(jj,2))~=0);
                            else
                                Temp = find(Fusion2(:,ChargeLabelsExternalUsed(jj,1),InverseIrrep(ChargeLabelsExternalUsed(jj,2)))~=0);
                            end
                            
                        else
                            Temp = find(Fusion2(:,InverseIrrep(ChargeLabelsExternalUsed(jj,1)),InverseIrrep(ChargeLabelsExternalUsed(jj,2)))~=0);
                        end
                        Temp = Temp(:);
                        if ~isempty(Temp)
                            
                            if numel(Temp)>1
                                ChargeLabelsExternalUsed(((jj+1):end)+numel(Temp)-1,:) = ChargeLabelsExternalUsed((jj+1):end,:);
                                ChargeLabelsExternalUsed((jj-1)+(1:numel(Temp)),:) = repmat(ChargeLabelsExternalUsed(jj,:),[numel(Temp),1]);
                                ChargeLabelsInternalUsed(((jj+1):end)+numel(Temp)-1,:) = ChargeLabelsInternalUsed((jj+1):end,:);
                                ChargeLabelsInternalUsed((jj-1)+(1:numel(Temp)),:) = repmat(ChargeLabelsInternalUsed(jj,:),[numel(Temp),1]);
                                for ll = 1:size(Temp,1);
                                    ChargeLabelsInternalUsed(jj+(ll-1)-1 + 1,1) = Temp(ll);
                                end
                            else
                                ChargeLabelsInternalUsed(jj,1) = Temp;
                            end
                        else
                            ChargeLabelsInternalUsed(jj,:) = [];
                            ChargeLabelsExternalUsed(jj,:) = [];
                        end
                        
                    end
                    
                    for kk = 1:(size(ChargeLabelsExternalUsed,2)-2)
                        
                        for jj = size(ChargeLabelsExternalUsed,1):-1:1
                            if ChargeDirectionsUsed(kk+2) ==1
                                Temp = find(Fusion2(:,ChargeLabelsInternalUsed(jj,kk),ChargeLabelsExternalUsed(jj,kk+2))~=0);
                            else
                                Temp = find(Fusion2(:,ChargeLabelsInternalUsed(jj,kk),InverseIrrep(ChargeLabelsExternalUsed(jj,kk+2)))~=0);
                            end
                            Temp = Temp(:);
                            if ~isempty(Temp)
                                
                                if numel(Temp)>1
                                    ChargeLabelsExternalUsed(((jj+1):end)+numel(Temp)-1,:) = ChargeLabelsExternalUsed((jj+1):end,:);
                                    ChargeLabelsExternalUsed((jj-1)+(1:numel(Temp)),:) = repmat(ChargeLabelsExternalUsed(jj,:),[numel(Temp),1]);
                                    ChargeLabelsInternalUsed(((jj+1):end)+numel(Temp)-1,:) = ChargeLabelsInternalUsed((jj+1):end,:);
                                    ChargeLabelsInternalUsed((jj-1)+(1:numel(Temp)),:) = repmat(ChargeLabelsInternalUsed(jj,:),[numel(Temp),1]);
                                    for ll = 1:size(Temp,1);
                                        ChargeLabelsInternalUsed(jj+(ll-1)-1 + 1,kk+1) = Temp(ll);
                                    end
                                else
                                    ChargeLabelsInternalUsed(jj,kk+1) = Temp;
                                end
                            else
                                ChargeLabelsInternalUsed(jj,:) = [];
                                ChargeLabelsExternalUsed(jj,:) = [];
                            end
                            
                        end
                        
                    end
                    
                    MultiplicitiesInternalUsed = ones([size(ChargeLabelsInternalUsed,1), TotalLegs-2]);
                
                    
                else %if multiplcitites are important
                    
                   %first one is special:
                    for jj = size(ChargeLabelsExternalUsed,1):-1:1
                        if ~any(ChargeDirectionsUsed([1,2]) ==-1)
                            
                            Temp = find(Fusion2(:,ChargeLabelsExternalUsed(jj,1),ChargeLabelsExternalUsed(jj,2))~=0);
                            
                        elseif ChargeDirectionsUsed(1)*ChargeDirectionsUsed(2) == -1
                            
                            if ChargeDirectionsUsed(1) == -1
                                Temp = find(Fusion2(:,InverseIrrep(ChargeLabelsExternalUsed(jj,1)),ChargeLabelsExternalUsed(jj,2))~=0);
                            else
                                Temp = find(Fusion2(:,ChargeLabelsExternalUsed(jj,1),InverseIrrep(ChargeLabelsExternalUsed(jj,2)))~=0);
                            end
                            
                        else
                            Temp = find(Fusion2(:,InverseIrrep(ChargeLabelsExternalUsed(jj,1)),InverseIrrep(ChargeLabelsExternalUsed(jj,2)))~=0);
                        end
                        Temp = Temp(:);
                        if ~isempty(Temp)
                            Copies = Fusion2(Temp,ChargeLabelsExternalUsed(jj,1),ChargeLabelsExternalUsed(jj,2));
                            if sum(Copies)>1
                                ChargeLabelsExternalUsed(((jj+1):end)+sum(Copies)-1,:) = ChargeLabelsExternalUsed((jj+1):end,:);
                                ChargeLabelsExternalUsed((jj-1)+(1:sum(Copies)),:) = repmat(ChargeLabelsExternalUsed(jj,:),[sum(Copies),1]);
                                ChargeLabelsInternalUsed(((jj+1):end)+sum(Copies)-1,:) = ChargeLabelsInternalUsed((jj+1):end,:);
                                ChargeLabelsInternalUsed((jj-1)+(1:sum(Copies)),:) = repmat(ChargeLabelsInternalUsed(jj,:),[sum(Copies),1]);
                                MultiplicitiesInternalUsed(((jj+1):end)+sum(Copies)-1,:) = MultiplicitiesInternalUsed((jj+1):end,:);
                                MultiplicitiesInternalUsed((jj-1)+(1:sum(Copies)),:) = repmat(MultiplicitiesInternalUsed(jj,:),[sum(Copies),1]);
                                for ll = 1:size(Temp,1);
                                    ChargeLabelsInternalUsed(jj+sum(Copies(1:(ll-1)))-1 + (1:Copies(ll)),1) = repmat(Temp(ll),[Copies(ll),1]);
                                    MultiplicitiesInternalUsed(jj+sum(Copies(1:(ll-1)))-1 + (1:Copies(ll)),1) = reshape(1:Copies(ll),[Copies(ll),1]);
                                end
                            else
                                ChargeLabelsInternalUsed(jj,1) = Temp;
                                MultiplicitiesInternalUsed(jj,1) = 1;
                            end
                        else
                            ChargeLabelsInternalUsed(jj,:) = [];
                            ChargeLabelsExternalUsed(jj,:) = [];
                        end
                        
                    end
                    
                    for kk = 1:(size(ChargeLabelsExternalUsed,2)-2)
                        
                        for jj = size(ChargeLabelsExternalUsed,1):-1:1
                            if ChargeDirectionsUsed(kk+2) ==1
                                Temp = find(Fusion2(:,ChargeLabelsInternalUsed(jj,kk),ChargeLabelsExternalUsed(jj,kk+2))~=0);
                            else
                                Temp = find(Fusion2(:,ChargeLabelsInternalUsed(jj,kk),InverseIrrep(ChargeLabelsExternalUsed(jj,kk+2)))~=0);
                            end 
                            Temp = Temp(:);
                            if ~isempty(Temp)
                                Copies = Fusion2(Temp,ChargeLabelsExternalUsed(jj,1),ChargeLabelsExternalUsed(jj,2));
                                if sum(Copies)>1
                                    ChargeLabelsExternalUsed(((jj+1):end)+sum(Copies)-1,:) = ChargeLabelsExternalUsed((jj+1):end,:);
                                    ChargeLabelsExternalUsed((jj-1)+(1:sum(Copies)),:) = repmat(ChargeLabelsExternalUsed(jj,:),[sum(Copies),1]);
                                    ChargeLabelsInternalUsed(((jj+1):end)+sum(Copies)-1,:) = ChargeLabelsInternalUsed((jj+1):end,:);
                                    ChargeLabelsInternalUsed((jj)+(1:sum(Copies)),:) = repmat(ChargeLabelsInternalUsed(jj,:),[sum(Copies),1]);
                                    MultiplicitiesInternalUsed(((jj+1):end)+sum(Copies)-1,:) = MultiplicitiesInternalUsed((jj+1):end,:);
                                    MultiplicitiesInternalUsed((jj)+(1:sum(Copies)),:) = repmat(MultiplicitiesInternalUsed(jj,:),[sum(Copies),1]);
                                    for ll = 1:size(Temp,1);
                                        ChargeLabelsInternalUsed(jj+sum(Copies(1:(ll-1)))-1 + (1:Copies(ll)),kk+1) = repmat(Temp(ll),[Copies(ll),1]);
                                        MultiplicitiesInternalUsed(jj+sum(Copies(1:(ll-1)))-1 + (1:Copies(ll)),kk+1) = reshape(1:Copies(ll),[Copies(ll),1]);
                                    end
                                else
                                    ChargeLabelsInternalUsed(jj,kk+1) = Temp;
                                    MultiplicitiesInternalUsed(jj,kk+1) = 1;
                                end
                            else
                                ChargeLabelsInternalUsed(jj,:) = [];
                                ChargeLabelsExternalUsed(jj,:) = [];
                                end
                                
                            end
                            
                        end
                        
                        MultiplicitiesInternalUsed = [MultiplicitiesInternalUsed];
                        
                end
                    
                end
            
            else
                
                ChargeLabelsInternalUsed = ones([size(ChargeLabelsExternalUsed,1),0]);
                MultiplicitiesInternalUsed = ones([size(ChargeLabelsExternalUsed,1),0]);
                
            end
        end
        
    end
    
    
    %these are the structure methods:
    
    methods(Access = 'public')
        
        function self = swap(self,j,k)
            
            AccessStructure = self.Structure;
            AccessChargeDirections = self.ChargeDirections;
            AccessStoredLocations = self.StoredLocations;
            AccessMultiplicitiesInternal = self.MultiplicitiesInternal;
            AccessChargeLabelsInternal = self.ChargeLabelsInternal;
            
            %check for errors:
            if (j>size(AccessStructure,2)) || (k>size(AccessStructure,2)) ||...
                    (j<0) || (k<0)
                error('These two don''t work for swapping, they aren''t in the correct reigme')
            end
            
            if (j~=round(j)) || (k~=round(k))
                error('one of these two aren''t integers')
            end
            
            
            %structure:
            StuctureTemp = AccessStructure(:,k+1);
            AccessStructure(:,k+1) = AccessStructure(:,j+1);
            StructureTemp(:,j+1) = StructureTemp;
            
            Locj = AccessStructure == j;
            Lock = AccessStructure == k;
            AccessStructure([Locj,Lock]) = [j,k];

            
            %ChargeDirections:
            ChargeTemp = AccessChargeDirections(:,k+1);
            AccessChargeDirections(:,k+1) = AccessChargeDirections(:,j+1);
            AccessChargeDirections(:,j+1) = ChargeTemp;
            
            %StoredLocations
            StoredLocationsTemp = AccessStoredLocations{j};
            AccessStoredLocations{j} = AccessStoredLocations{k};
            AccessStoredLocations{k} = StoredLocationsTemp;
            
            
            %ChargeLabelsInternal
            ChargeLabelsInternalTemp = AccessChargeLabelsInternal(:,j);
            AccessChargeLabelsInternal(:,j) = AccessChargeLabelsInternal(:,k);
            AccessChargeLabelsInternal(:,k) = ChargeLabelsInternalTemp;
            
            
            %MultiplicitiesInternal
            MultiplicitiesInternalTemp = AccessMultiplicitiesInternal(:,j);
            AccessMultiplicitiesInternal(:,j) = AccessMultiplicitiesInternal(:,k);
            AccessMultiplicitiesInternal(:,k) = MultiplicitiesInternalTemp;
            
            
            self.Structure = AccessStructure;
            self.ChargeDirections = AccessChargeDirections;
            self.StoredLocations = AccessStoredLocations;
            self.MultiplicitiesInternal = AccessMultiplicitiesInternal;
            self.ChargeLabelsInternal = AccessChargeLabelsInternal;
            
        end
        
        function self = FMove(self,Position,MoveDirection)
            
            
            %Things we need to change:
            % - Structure
            % - ChargeDirections
            % - StoredLocations
            %
            % - ChargeLabelsInternal
            % - MultiplicitiesInternal
            % - TensorEntries
            
            
            AccessStructure = self.Structure;
            AccessChargeDirections = self.ChargeDirections;
            AccessStoredLocations = self.StoredLocations;
            
            AccessMultiplicitiesInternal = self.MultiplicitiesInternal;
            AccessChargeLabelsInternal = self.ChargeLabelsInternal;
            AccessChargeLabelsExternal = self.ChargeLabelsExternal;
            AccessTensorEntries = self.TensorEntries;
            
            %checks:
            
            if (Position>size(AccessStructure,2)) || (Position<0)
                error('This Position value doesn''t work for F-moves, This should be the top label')
            end
            
            if (Position~=round(Position))
                error('The Position given for the F-move isn''t an integer')
            end
            
            if ~((MoveDirection==1)||(MoveDirection==-1))
                error('The MoveDirection parameter should be +1 (for moving the leg left to right) or -1 (for moving the leg right to left)')
            end
            
            FlagMovingToRight = MoveDirection == 1;
            
            if AccessStructure(2-FlagMovingToRight,Position+1)<0
                error('This is not a valid Fmove, you''re trying to use an external leg for the middle leg')
            end
            
            
            %first do the degeneracy/charge independent stuff
            
            if FlagMovingToRight
                
                N = AccessStructure(1,Position+1);
                C = AccessStructure(2,Position+1);
                A = AccessStructure(1,N+1);
                B = AccessStructure(2,N+1);
                
                AccessStructure(1,Position+1) = A;
                AccessStructure(2,Position+1) = N;
                AccessStructure(1,N+1) = B;
                AccessStructure(2,N+1) = C;
                
                CDN = AccessChargeDirections(1,Position+1);
                CDC = AccessChargeDirections(2,Position+1);
                CDA = AccessChargeDirections(1,N+1);
                CDB = AccessChargeDirections(2,N+1);
                
                AccessChargeDirections(1,Position+1) = CDA;
                AccessChargeDirections(2,Position+1) = CDN;
                AccessChargeDirections(1,N+1) = CDB;
                AccessChargeDirections(2,N+1) = CDC;
                
                if B>0
                    StoredLeft = AccessStoredLocations{B};
                else
                    StoredLeft = -B;
                end
                
                if C>0
                    StoredRight = AccessStoredLocations{C};
                else
                    StoredRight = -C;
                end
                
                AccessStoredLocations{N} = [StoredLeft,StoredRight];
                
            else
                
                A = AccessStructure(1,Position+1);
                N = AccessStructure(2,Position+1);
                B = AccessStructure(1,N+1);
                C = AccessStructure(2,N+1);
                
                AccessStructure(1,Position+1) = N;
                AccessStructure(2,Position+1) = C;
                AccessStructure(1,N+1) = A;
                AccessStructure(2,N+1) = B;
                
                CDA = AccessChargeDirections(1,Position+1);
                CDN = AccessChargeDirections(2,Position+1);
                CDB = AccessChargeDirections(1,N+1);
                CDC = AccessChargeDirections(2,N+1);
                
                AccessChargeDirections(1,Position+1) = CDN;
                AccessChargeDirections(2,Position+1) = CDC;
                AccessChargeDirections(1,N+1) = CDA;
                AccessChargeDirections(2,N+1) = CDB;
                
                if A>0
                    StoredLeft = AccessStoredLocations{A};
                else
                    StoredLeft = -A;
                end
                
                if B>0
                    StoredRight = AccessStoredLocations{B};
                else
                    StoredRight = -B;
                end
                
                AccessStoredLocations{N} = [StoredLeft,StoredRight];
                
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now for applying Fs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            LocA = abs(A)+(size(AccessChargeLabelsExternal,2))*(A>0);
            LocB = abs(B)+(size(AccessChargeLabelsExternal,2))*(B>0);
            LocC = abs(C)+(size(AccessChargeLabelsExternal,2))*(C>0);
            
            LocN = N+size(AccessChargeLabelsExternal,2);%has to be an internal by definition.
            LocD = Position+size(AccessChargeLabelsExternal,2);
            
            ChargesTotal = [AccessChargeLabelsExternal, AccessChargeLabelsInternal];
            
            BlockLocations = [LocA,LocB,LocC,LocD];
            ListLocations = [LocN];
            
            if self.SymHandle.IsMultiplicities
                LocMu = size(AccessChargeLabelsExternal,2)+size(AccessChargeLabelsInternal,2)+N;
                LocNu = size(AccessChargeLabelsExternal,2)+size(AccessChargeLabelsInternal,2)+Position;
                ChargesTotal = [ChargesTotal, AccessMultiplicitiesInternal];
                ListLocations = [LocN,LocMu, LocNu];
            end
            
            ListLocationsComplimentary = 1:size(ChargesTotal,2);
            ListLocationsComplimentary([ListLocations,BlockLocations]) =[];
            
           
            
            [UniqueBlocks,~,IndexFBlocks] = unique(ChargesTotal(:,BlockLocations),'rows');
                
            [F,CumSumF,FLabels] = self.SymHandle.getFMoves;
                   
            
            if ~self.SymHandle.IsNonAbelian 
                
                %Abelian case
                %make a cell of size of the unique blocks
                
                %Just need to multiply the matrix
                
                Numbers = 1:numel(IndexFBlocks);
                
                for jj = 1:size(UniqueBlocks,1)
                    
                    %here we are fixing a b c d as they won't mix due to the
                    %block diagonal structure of F matrix:
                    
                    FLocation = sum((UniqueBlocks(jj,:)-1).*(size(self.SymHandle.getDim,1).^[0,1,2]))+1;
                    
                    FBlock = F{FLocation,2-FlagMovingToRight};
                    
                    Keep = Numbers(IndexFBlocks == jj);
                    
                    %Multiply by relevent F Note that because we have fixed
                    %the blocks we have also fixed the Old and new labels
                    %(Because this is abelian):
                    if FBlock ~= 1 %Special Simple case
                        for kk = Keep
                            AccessTensorEntries{kk} = AccessTensorEntries{kk}*FBlock;
                        end
                    end
                    
                    AccessChargeLabelsInternal(N,Keep) = repmat(FLabels{FLocation, 1+FlagMovingToRight},[1,numel(Keep)]);
                end
                
                
            else
                %if it is non-abelian
                
                if ~self.SymHandle.IsMultiplicities 
                    %if there are no multiplicities
                    
                    %make a cell of size of the unique blocks
                    
                    TensorEntriesNew = cell([size(UniqueBlocks,1),1]);
                    ChargesNew = cell([size(UniqueBlocks,1),1]);
                    
                    
                    for jj = 1:size(UniqueBlocks,1)
                        
                        %here we are fixing a b c d as they won't mix due to the
                        %block diagonal structure of F matrix:
                        
                        
                        Keep = (IndexFBlocks == jj);
                        
                        TensorEntriesBlock = AccessTensorEntries(Keep);
                                            
                        [UniqueFCaring,~,IndexFCaring] = unique(ChargesTotal(Keep,ListLocations),'rows');
                        [UniqueFUnCaring,FirstFUnCaring,IndexFUnCaring] = unique(ChargesTotal(Keep,ListLocationsComplimentary),'rows');
                        Dimensions = ones([size(FirstFUnCaring,1),size(AccessChargeLabelsExternal,2)]);
                        
                        %note that we have a minimum of 3 dimensions so we
                        %don't need to consider special cases of 0 or 1
                        %legs
                        for kk = 1:size(FirstFUnCaring,1)
                            Temp = size(TensorEntriesBlock{kk});
                            Dimensions(kk,1:length(Temp)) = Temp;
                        end
                        DimensionsSingle = prod(Dimensions,2)';
                        
                        for kk = 1:numel(TensorEntriesBlock)
                            TensorEntriesBlock{kk} = reshape(TensorEntriesBlock{kk}, [1,numel(AccessTensorEntries{kk})]);
                        end
                        
                        TensorEntriesBlock2 = cell([size(UniqueFCaring,1),size(UniqueFUnCaring,1)]);
                        
                        TensorEntriesBlock2(IndexFCaring+(IndexFUnCaring-1)*size(IndexFCaring,1)) = TensorEntriesBlock(:);
                        
                        %do the thing like the function to see if a cell is
                        %empty.
                        
                        EmptyEntries = find(  (cellfun('isempty', TensorEntriesBlock2(:))));
                        
                        kk1 = mod(EmptyEntries-1,size(TensorEntriesBlock2,1))+1;
                        kk2 = round((EmptyEntries-kk1)/size(TensorEntriesBlock2,1));
                        kk1 = kk1(:); kk2 = kk2(:)+1;
                        
                        for kk = [kk1,kk2]'
                            TensorEntriesBlock2{kk(1),kk(2)} = zeros([1,DimensionsSingle(kk(2))]);
                        end
                        TensorEntriesBlock = cell2mat(TensorEntriesBlock2);
                        
                        %multiply by relevent F block
                        
                        
                        FLocation = sum((UniqueBlocks(jj,:)-1).*(size(self.SymHandle.getDim,1).^[0,1,2,3]))+1;
                        FBlock = F{FLocation,2-FlagMovingToRight};
                        
                           
                        LocationInput = zeros([size(UniqueFCaring,1),1]);
                        for kk = 1:size(UniqueFCaring,1)
                            LocationInput(kk) = find(FLabels{FLocation,2-FlagMovingToRight}==UniqueFCaring(kk)); 
                            %recall that we only have a single entry in
                            %UniqueFCaring as we don't worry about
                            %multiplicities and only care about one leg
                        
                        end
                        FBlock = FBlock(:,LocationInput);
                        
                        %now work out what the outputs have to be (which
                        %has been pre-computed for us).
                        OutputRelabeled = FLabels{FLocation,1+FlagMovingToRight};
                        
                        Temp = any(FBlock~=0,2)';
                        OutputRelabeled = OutputRelabeled(Temp,:);
                        FBlock = FBlock(Temp,:);
                            
                        
                        TensorEntriesBlock = FBlock*TensorEntriesBlock;
                        
                        %pull apart using DimensionsSingle
                        
                        UsefulBlocks = [repmat(1:size(TensorEntriesBlock,1),[1,numel(DimensionsSingle)]);...
                            reshape(repmat((1:numel(DimensionsSingle))',[1,size(TensorEntriesBlock,1)]),[1,numel(DimensionsSingle)*size(TensorEntriesBlock,1)])];
    
                        %for the moment we will assume they are all filled so
                        %we won't bother working out which ones are empty.
                        
                        TensorEntriesBlock = mat2cell(TensorEntriesBlock,ones([1,size(TensorEntriesBlock,1)]),DimensionsSingle);
                        
                        for kk = UsefulBlocks
                            TensorEntriesBlock{kk(1),kk(2)} = reshape(TensorEntriesBlock{kk(1),kk(2)}, Dimensions(kk(2),:));
                        end
                        
                        
                        %add to the end of the tensor list.
                        TensorEntriesNew{jj} = TensorEntriesBlock(:);
                        
                        
                        ChargesNew{jj} = [repmat(UniqueBlocks(jj,:),[numel(TensorEntriesBlock),1]),...
                            OutputRelabeled,...
                            UniqueFUnCaring(reshape(repmat(1:size(TensorEntriesBlock,2),[size(TensorEntriesBlock,1),1]),[1,numel(TensorEntriesBlock)]),:)];
                        
                    end
                     
                    
                else
                    
                    
                    
                    %make a cell of size of the unique blocks
                    
                    TensorEntriesNew = cell([size(UniqueBlocks,1),1]);
                    ChargesNew = cell([size(UniqueBlocks,1),1]);
                    
                    
                    for jj = 1:size(UniqueBlocks,1)
                        
                        %here we are fixing a b c d as they won't mix due to the
                        %block diagonal structure of F matrix:
                        
                        
                        Keep = (IndexFBlocks == jj);
                        
                        TensorEntriesBlock = AccessTensorEntries(Keep);
                                            
                        [UniqueFCaring,~,IndexFCaring] = unique(ChargesTotal(Keep,ListLocations),'rows');
                        [UniqueFUnCaring,FirstFUnCaring,IndexFUnCaring] = unique(ChargesTotal(Keep,ListLocationsComplimentary),'rows');
                        Dimensions = ones([size(FirstFUnCaring,1),size(AccessChargeLabelsExternal,2)]);
                        for kk = 1:size(FirstFUnCaring,1)
                            Temp = size(TensorEntriesBlock{kk});
                            Dimensions(kk,1:length(Temp)) = Temp;
                        end
                        DimensionsSingle = prod(Dimensions,2)';
                        
                        for kk = 1:numel(TensorEntriesBlock)
                            TensorEntriesBlock{kk} = reshape(TensorEntriesBlock{kk}, [1,numel(AccessTensorEntries{kk})]);
                        end
                        
                        TensorEntriesBlock2 = cell([size(UniqueFCaring,1),size(UniqueFUnCaring,1)]);
                        
                        TensorEntriesBlock2(IndexFCaring+(IndexFUnCaring-1)*size(IndexFCaring,1)) = TensorEntriesBlock(:);
                        
                        %do the thing like the function to see if a cell is
                        %empty.
                        
                        EmptyEntries = find(  (cellfun('isempty', TensorEntriesBlock2(:))));
                        
                        kk1 = mod(EmptyEntries-1,size(TensorEntriesBlock2,1))+1;
                        kk2 = round((EmptyEntries-kk1)/size(TensorEntriesBlock2,1));
                        kk1 = kk1(:); kk2 = kk2(:)+1;
                        
                        for kk = [kk1,kk2]'
                            TensorEntriesBlock2{kk(1),kk(2)} = zeros([1,DimensionsSingle(kk(2))]);
                        end
                        
                        
                        %multiply by relevent F block
                        
                        
                        FLocation = sum((UniqueBlocks(jj,:)-1).*(size(self.SymHandle.getDim,1).^[0,1,2,3]))+1;
                        FBlock = F{FLocation,2-FlagMovingToRight};
                        
                        
                            
                        
                        if ~self.SymHandle.IsMultiplicities
                            %if we don't need to worry about multiplicities
                            
                            %first need to convert to 
                            
                            LocationInput = zeros([size(UniqueFCaring,1),1]);
                            for kk = 1:size(UniqueFCaring,1)
                                LocationInput(kk) = find(FLabels{FLocation,2-FlagMovingToRight}==UniqueFCaring(kk)); 
                                %recall that we only have a single entry in
                                %UniqueFCaring as we don't worry about
                                %multiplicities and only care about one leg
                                
                            end
                            
                            FBlock = FBlock(:,LocationInput);
                            
                            %now work out what the outputs have to be (which
                            %has been pre-computed for us).
                            OutputRelabeled = FLabels{FLocation,1+FlagMovingToRight};
                            
                            
                            Temp = any(FBlock~=0,2)';
                            OutputRelabeled = OutputRelabeled(Temp,:);
                            FBlock = FBlock(Temp,:);
                            
                        else
                            %if we need to worry about multiplicities
                            
                            %STILL DOING AT THE MOMENT
                            
                        end
                    
                    
                        %pull apart using DimensionsSingle
                        
                        UsefulBlocks = [repmat(1:size(TensorEntriesBlock,1),[1,numel(DimensionsSingle)]);...
                            reshape(repmat((1:numel(DimensionsSingle))',[1,size(TensorEntriesBlock,1)]),[1,numel(DimensionsSingle)*size(TensorEntriesBlock,1)])];
    
                        %for the moment we will assume they are all filled so
                        %we won't bother working out which ones are empty.
                        
                        TensorEntriesBlock = mat2cell(TensorEntriesBlock,ones([1,size(TensorEntriesBlock,1)]),DimensionsSingle);
                        
                        
                        
                        for kk = UsefulBlocks
                            TensorEntriesBlock{kk(1),kk(2)} = reshape(TensorEntriesBlock{kk(1),kk(2)}, Dimensions(kk(2),:));
                        end
                        
                        
                        %add to the end of the tensor list.
                        TensorEntriesNew{jj} = TensorEntriesBlock(:);
                        
                        
                        ChargesNew{jj} = [repmat(UniqueBlocks(jj,:),[numel(TensorEntriesBlock),1]),...
                            OutputRelabeled,...
                            UniqueFUnCaring(reshape(repmat(1:size(TensorEntriesBlock,2),[size(TensorEntriesBlock,1),1]),[1,numel(TensorEntriesBlock)]),:)];
                        
                    end
                    
                end
                
                
                AccessTensorEntries = vertcat(TensorEntriesNew{:});
                ChargesNew = cell2mat(ChargesNew);
                                
                %reorganise
                [~,Index] = sort([BlockLocations,ListLocations,ListLocationsComplimentary]);
                ChargesNew = ChargesNew(:,Index);
                
                
                AccessChargeLabelsExternal = ChargesNew(:,1:size(AccessChargeLabelsExternal,2));
                AccessChargeLabelsInternal = ChargesNew(:,size(AccessChargeLabelsExternal,2)+(1:size(AccessChargeLabelsInternal,2)));
                
                if self.SymHandle.IsMultiplicities
                    AccessMultiplicitiesInternal = ChargesNew(:,(size(AccessChargeLabelsExternal)+size(AccessChargeLabelsInternal,2)+1):end);
                end
                %combine the tensor lists together (doing here so we don't
                %have to worry about the cost of extending the list)
                    
                    
            end
            
            
            
            
            
            
            self.Structure = AccessStructure;
            self.ChargeDirections = AccessChargeDirections;
            self.StoredLocations = AccessStoredLocations;
            
            self.MultiplicitiesInternal = AccessMultiplicitiesInternal;
            self.ChargeLabelsInternal = AccessChargeLabelsInternal;
            self.ChargeLabelsExternal = AccessChargeLabelsExternal;
            self.TensorEntries = AccessTensorEntries;
            
        end
        
        function self = RMove(self,Position, MoveDirection)
            %Things we need to change:
            % - Structure
            % - ChargeDirections
            % - Braidings
            % - BraidingDirections
            %
            % - MultiplicitiesInternal
            % - TensorEntries
            
            
            AccessStructure = self.Structure;
            AccessChargeDirections = self.ChargeDirections;
            AccessBraidings = self.Braidings;
            AccessBraidingDirections = self.BraidingDirections;
            
            AccessMultiplicitiesInternal = self.MultiplicitiesInternal;
            AccessChargeLabelsInternal = self.ChargeLabelsInternal;
            AccessChargeLabelsExternal = self.ChargeLabelsExternal;
            AccessTensorEntries = self.TensorEntries;
            
            
            AccessStoredLocations = self.StoredLocations;
            
            %checks:
            
            if (Position>size(AccessStructure,2)) || (Position<=0)
                error('This Position value doesn''t work for R-moves, This should be the c label')
            end
            
            if (Position~=round(Position))
                error('The Position given for the R-move isn''t an integer')
            end
            
            if ~((MoveDirection==1)||(MoveDirection==-1))
                error('The MoveDirection parameter should be +1 (for left leg over right leg) or -1 (for right leg over left leg)')
            end
            
            FlagLeftOverRight = MoveDirection == 1;
            
            
            %first do the degeneracy/charge independent stuff
            
            AccessStructure([2,1], Position+1) = AccessStructure(:,Position+1);
            AccessChargeDirections([2,1], Position+1) = AccessChargeDirections(:,Position+1);
            
            
            
            %the swaps always happen with the one going over happening
            %first one by one
            
            if AccessStructure(1,Position+1)<0
                LeftStoredLocations = -AccessStructure(1,Position+1);
                LocLeft = -AccessStructure(1,Position+1);
            else
                LeftStoredLocations = AccessStoredLocations{AccessStructure(1,Position+1)};
                LocLeft = AccessStructure(1,Position+1) + size(AccessChargeLabelsExternal,2);
            end
            
            if AccessStructure(2,Position+1)<0
                RightStoredLocations = -AccessStructure(2,Position+1);
                LocRight = -AccessStructure(2,Position+1);
            else
                RightStoredLocations = AccessStoredLocations{AccessStructure(2,Position+1)};
                LocRight = AccessStructure(2,Position+1) + size(AccessChargeLabelsExternal,2);
            end
            
            
            if FlagLeftOverRight
                
                Swaps = [reshape(repmat(LeftStoredLocations,[size(RightStoredLocations,2),1]),[1,size(RightStoredLocations,2)*size(LeftStoredLocations,2)]);...
                    repmat(RightStoredLocations,[1,size(LeftStoredLocations,2)])];
                
            else
                
                Swaps = [repmat(LeftStoredLocations,[1,size(RightStoredLocations,2)]);...
                    reshape(repmat(RightStoredLocations,[size(LeftStoredLocations,2),1]),[1,size(RightStoredLocations,2)*size(LeftStoredLocations,2)])];
                
            end
            
            AccessBraidings = [AccessBraidings,Swaps];
            AccessBraidingDirections = [AccessBraidingDirections, (1-2*FlagLeftOverRight)*ones([1,size(Swaps,2)])];
            
            
            
            
            
            
            
            
            %LocA = abs(A)+(size(AccessChargeLabelsExternal,2))*(A>0);
            %LocB = abs(B)+(size(AccessChargeLabelsExternal,2))*(B>0);
            %LocC = abs(C)+(size(AccessChargeLabelsExternal,2))*(C>0);
            
            ChargesTotal = [AccessChargeLabelsExternal, AccessChargeLabelsInternal];
            
            BlockLocations = [LocLeft,LocRight,Position+size(AccessChargeLabelsExternal,2)]; %[LocA,LocB,LocC];
            
            BlockLocationsOut = BlockLocations([2,1,3]);
            
            ListLocations = zeros([1,0]);
            
            if self.SymHandle.IsMultiplicities
                ListLocations = size(AccessChargeLabelsExternal,2)+size(AccessChargeLabelsInternal,2)+Position;
            end
            
            
            ListLocationsComplimentary = 1:size(ChargesTotal,2);
            ListLocationsComplimentary([ListLocations,BlockLocations]) =[];
            
           
            
            R = self.SymHandle.getRMoves;%unlike in the F move case our blocks are square with the multiplicity as the numbe of degrees of freedom.
            
            
            if ~self.SymHandle.IsMultiplicities %If no muliplicities, either abelian or non-abelian.
                
                [UniqueBlocks,~,IndexFBlocks] = unique(ChargesTotal(:,BlockLocations),'rows');
                
                
                %Just need to multiply the matrix
                
                Numbers = 1:numel(IndexFBlocks);
                
                for jj = 1:size(UniqueBlocks,1)
                    
                    %here we are fixing a b c as they won't mix due to the
                    %block diagonal structure of R matrix:
                    
                    RLocation = sum((UniqueBlocks(jj,:)-1).*(size(self.SymHandle.getDim,1).^[0,1,2]))+1;
                    
                    RBlock = R{RLocation,2-FlagLeftOverRight};
                    
                    %Multiply by relevent R
                    if RBlock ~=1
                        Keep = Numbers(IndexFBlocks == jj);
                        
                        for kk = Keep
                            AccessTensorEntries{kk} = AccessTensorEntries{kk}*RBlock;
                        end
                    end
                    
                end
                 
                
                
            else
                
                
                %make a cell of size of the unique blocks
                
                TensorEntriesNew = cell([size(UniqueBlocks,1),1]);
                ChargesNew = cell([size(UniqueBlocks,1),1]);
                
                
                for jj = 1:size(UniqueBlocks,1)
                    
                    %here we are fixing a b c d as they won't mix due to the
                    %block diagonal structure of F matrix:
                    
                    
                    Keep = (IndexFBlocks == jj);
                    
                    TensorEntriesBlock = AccessTensorEntries(Keep);
                                        
                    [UniqueFCaring,~,IndexFCaring] = unique(ChargesTotal(Keep,ListLocations),'rows');
                    [UniqueFUnCaring,FirstFUnCaring,IndexFUnCaring] = unique(ChargesTotal(Keep,ListLocationsComplimentary),'rows');
                    Dimensions = ones([size(FirstFUnCaring,1),size(AccessChargeLabelsExternal,2)]);
                    for kk = 1:size(FirstFUnCaring,1)
                        Temp = size(TensorEntriesBlock{kk});
                        Dimensions(kk,1:length(Temp)) = Temp;
                    end
                    DimensionsSingle = prod(Dimensions,2)';
                    
                    for kk = 1:numel(TensorEntriesBlock)
                        TensorEntriesBlock{kk} = reshape(TensorEntriesBlock{kk}, [1,numel(AccessTensorEntries{kk})]);
                    end
                    
                    TensorEntriesBlock2 = cell([size(UniqueFCaring,1),size(UniqueFUnCaring,1)]);
                    
                    TensorEntriesBlock2(IndexFCaring+(IndexFUnCaring-1)*size(IndexFCaring,1)) = TensorEntriesBlock(:);
                    
                    %do the thing like the function to see if a cell is
                    %empty.
                    
                    EmptyEntries = find(  (cellfun('isempty', TensorEntriesBlock2(:))));
                    
                    kk1 = mod(EmptyEntries-1,size(TensorEntriesBlock2,1))+1;
                    kk2 = round((EmptyEntries-kk1)/size(TensorEntriesBlock2,1));
                    kk1 = kk1(:); kk2 = kk2(:)+1;
                    
                    for kk = [kk1,kk2]'
                        TensorEntriesBlock2{kk(1),kk(2)} = zeros([1,DimensionsSingle(kk(2))]);
                    end
                    
                    
                    %multiply by relevent F block
                    
                    
                    FLocation = sum((UniqueBlocks(jj,:)-1).*(size(self.SymHandle.getDim,1).^[0,1,2,3]))+1;
                    FBlock = F{FLocation,2-FlagMovingToRight};
                    
                    
                        
                    
                    if ~self.SymHandle.IsMultiplicities
                        %if we don't need to worry about multiplicities
                        
                        %first need to convert to 
                        
                        LocationInput = zeros([size(UniqueFCaring,1),1]);
                        for kk = 1:size(UniqueFCaring,1)
                            LocationInput(kk) = find(FLabels{FLocation,2-FlagMovingToRight}==UniqueFCaring(kk)); 
                            %recall that we only have a single entry in
                            %UniqueFCaring as we don't worry about
                            %multiplicities and only care about one leg
                            
                        end
                        
                        FBlock = FBlock(:,LocationInput);
                        
                        %now work out what the outputs have to be (which
                        %has been pre-computed for us).
                        OutputRelabeled = FLabels{FLocation,1+FlagMovingToRight};
                        
                        
                        Temp = any(FBlock~=0,2)';
                        OutputRelabeled = OutputRelabeled(Temp,:);
                        FBlock = FBlock(Temp,:);
                        
                    else
                        %if we need to worry about multiplicities
                        
                        %STILL DOING AT THE MOMENT
                        
                    end
                    
                    
                    %pull apart using DimensionsSingle
                    
                    UsefulBlocks = [repmat(1:size(TensorEntriesBlock,1),[1,numel(DimensionsSingle)]);...
                        reshape(repmat((1:numel(DimensionsSingle))',[1,size(TensorEntriesBlock,1)]),[1,numel(DimensionsSingle)*size(TensorEntriesBlock,1)])];

                    %for the moment we will assume they are all filled so
                    %we won't bother working out which ones are empty.
                    
                    TensorEntriesBlock = mat2cell(TensorEntriesBlock,ones([1,size(TensorEntriesBlock,1)]),DimensionsSingle);
                    
                    
                    
                    for kk = UsefulBlocks
                        TensorEntriesBlock{kk(1),kk(2)} = reshape(TensorEntriesBlock{kk(1),kk(2)}, Dimensions(kk(2),:));
                    end
                    
                    
                    %add to the end of the tensor list.
                    TensorEntriesNew{jj} = TensorEntriesBlock(:);
                    
                    
                    ChargesNew{jj} = [repmat(UniqueBlocks(jj,:),[numel(TensorEntriesBlock),1]),...
                        OutputRelabeled,...
                        UniqueFUnCaring(reshape(repmat(1:size(TensorEntriesBlock,2),[size(TensorEntriesBlock,1),1]),[1,numel(TensorEntriesBlock)]),:)];
                    
                end
                 
                AccessTensorEntries = vertcat(TensorEntriesNew{:});
                AccessTensorEntries = vertcat(AccessTensorEntries{:});
                ChargesNew = cell2mat(ChargesNew);
                
                
                %reorganise
                [~,Index] = sort([BlockLocations,ListLocations,ListLocationsComplimentary]);
                ChargesNew = ChargesNew(:,Index);
                
                
                AccessChargeLabelsExternal = ChargesNew(:,1:size(AccessChargeLabelsExternal,2));
                AccessChargeLabelsInternal = ChargesNew(:,size(AccessChargeLabelsExternal,2)+(1:size(AccessChargeLabelsInternal,2)));
                
                if self.SymHandle.IsMultiplicities
                    AccessMultiplicitiesInternal = ChargesNew(:,(size(AccessChargeLabelsExternal)+size(AccessChargeLabelsInternal,2)+1):end);
                end
                %combine the tensor lists together (doing here so we don't
                %have to worry about the cost of extending the list)
                
                
                
                
                
            end
            
            
            
            
            
            
            
            self.Structure = AccessStructure;
            self.ChargeDirections = AccessChargeDirections;
            self.Braidings = AccessBraidings;
            self.BraidingDirections = AccessBraidingDirections;
            
            self.MultiplicitiesInternal = AccessMultiplicitiesInternal;
            self.TensorEntries = AccessTensorEntries;
            
        end
        
        function Bool = SymIsEqual(A,B, AllowedError)
            
            if ~isa(B,'SymTensor')
                error('A and B must both be SymTensors')
            end
            
            if ~isa(A,'SymTensor')
                error('A and B must both be SymTensors')
            end
            
            Bool = true;
            
            if nargin<3
                AllowedError = 10^-14;
            end
            
            A = A.SortLabels;
            B = B.SortLabels;
            
            if ~isequal(A.Structure,B.Structure)
                Bool = false;
                return
            end
            if ~isequal(A.ChargeDirections,B.ChargeDirections)
                Bool = false;
                return
            end
            
            ChargeLabelsA = [A.ChargeLabelsExternal,A.ChargeLabelsInternal,A.MultiplicitiesInternal];
            ChargeLabelsB = [B.ChargeLabelsExternal,B.ChargeLabelsInternal,B.MultiplicitiesInternal];
            
            if ~isequal(A.ChargeLabelsExternal,B.ChargeLabelsExternal)
                Bool = false;
                return
            end
            if ~isequal(A.ChargeLabelsInternal,B.ChargeLabelsInternal)
                Bool = false;
                return
            end
            if ~isequal(A.MultiplicitiesInternal,B.MultiplicitiesInternal)
                Bool = false;
                return
            end
            if ~isequal(A.Braidings,B.Braidings)
                Bool = false;
                return
            end
            if ~isequal(A.BraidingDirections,B.BraidingDirections)
                Bool = false;
                return
            end
            if ~isequal(A.ChargeLegDimensions,B.ChargeLegDimensions)
                Bool = false;
                return
            end
            if ~isequal(A.SymHandle,B.SymHandle)
                Bool = false;
                return
            end
            
            %if everything else is the same then the sizes of the tensor
            %parts and the components have to be the same
            for kk = 1:numel(A.TensorEntries)
                Diff = A.TensorEntries{kk} - B.TensorEntries{kk};
                if max(abs(Diff(:)))>AllowedError;
                    Bool = false;
                    return;
                end
            end
            
        end
        
        function Bool = SymIsSameShape(A,B)
            
            if ~isa(B,'SymTensor')
                error('A and B must both be SymTensors')
            end
            
            if ~isa(A,'SymTensor')
                error('A and B must both be SymTensors')
            end
            
            Bool = true;
            
            A = A.SortLabels;
            B = B.SortLabels;
            
            if ~isequal(A.Structure,B.Structure)
                Bool = false;
                return
            end
            if ~isequal(A.ChargeDirections,B.ChargeDirections)
                Bool = false;
                return
            end
            
            if ~isequal(A.Braidings,B.Braidings)
                Bool = false;
                return
            end
            if ~isequal(A.BraidingDirections,B.BraidingDirections)
                Bool = false;
                return
            end
            if ~isequal(A.SymHandle,B.SymHandle)
                Bool = false;
                return
            end
            
            if (~isequal(size(A.ChargeLegDimensions),size(B.ChargeLegDimensions))||~isequal(A.SymHandle,B.SymHandle))
                Bool = false;
                return
            elseif ~all(all(A.ChargeLegDimensions==B.ChargeLegDimensions|A.ChargeLegDimensions==0|B.ChargeLegDimensions==0))
                Bool = false;
                return
            end
        end
        
        function Bool = SymCanInsert(A,B)
            
            if ~isa(B,'SymTensor')
                error('A and B must both be SymTensors')
            end
            
            if ~isa(A,'SymTensor')
                error('A and B must both be SymTensors')
            end
            
            Bool = true;
            
            A = A.SortLabels;
            B = B.SortLabels;
            
            if ~isequal(A.Structure,B.Structure)
                Bool = false;
                return
            end
            if ~isequal(A.ChargeDirections,B.ChargeDirections)
                Bool = false;
                return
            end
            
            if ~isequal(A.Braidings,B.Braidings)
                Bool = false;
                return
            end
            if ~isequal(A.BraidingDirections,B.BraidingDirections)
                Bool = false;
                return
            end
            if ~isequal(A.SymHandle,B.SymHandle)
                Bool = false;
                return
            end
            
            if (~isequal(size(A.ChargeLegDimensions),size(B.ChargeLegDimensions))||~isequal(A.SymHandle,B.SymHandle))
                Bool = false;
                return
            elseif ~all(all(A.ChargeLegDimensions==B.ChargeLegDimensions|A.ChargeLegDimensions==0|B.ChargeLegDimensions==0))
                Bool = false;
                return
            end
            
            ChargeLabelsA = [A.ChargeLabelsExternal,A.ChargeLabelsInternal,A.MultiplicitiesInternal];
            ChargeLabelsB = [B.ChargeLabelsExternal,B.ChargeLabelsInternal,B.MultiplicitiesInternal];
            
            [ChargeLabelsUnique,~,~] = unique([ChargeLabelsB;ChargeLabelsA],'rows');
            
            if size(ChargeLabelsUnique,1) == size(ChargeLabelsB,1)
                %all good
            elseif size(ChargeLabelsUnique,1) > size(ChargeLabelsB,1)
                Bool = false;
                return
            else %size(ChargeLabelsUnique,1) < size(ChargeLabelsA,1)
                error('Error: Somehow Tensor B does not have unique ChargeLabels')
            end
            
        end
        
        function self = Permute(self,PermuteOrder)
            
            %HERAFix need to add checks
            
%            warning('Warning: Should use FreePermute, Permute currently isn''t coded to correct for internal sturcture yet');
            
            for kk = 1:numel(self.TensorEntries)
                self.TensorEntries{kk} = permute(self.TensorEntries{kk},PermuteOrder);
            end
        end
        
    end
    
end


function bool = IsInteger(A)
    bool = isnumeric(A);
    bool = bool&~any(abs(A(:)-round(A(:)))>10^-14);
end
