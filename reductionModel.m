function matPredict=reductionModel(Kd,Km,matOrig,Option,p)
         [Ud,ed]=eig(Kd);
         e=min(diag(ed));
         if e<0
            Kd=Kd-e*eye(size(Kd,1));
            [Ud,ed]=eig(Kd);
         end
         
         [Um,em]=eig(Km);
         e=min(diag(em));
         if e<0
            Km=Km-e*eye(size(Km,1));
            [Um,em]=eig(Km);
         end
         
         if Option==1
            z=reshape(Um'*matOrig'*Ud,numel(Um'*matOrig'*Ud),1);
         
            Ed=diag(ed);
            Ed=repmat(Ed',size(em,1),1);
            Ed=reshape(Ed,numel(Ed),1);
            Em=diag(em);
            Em=repmat(Em,size(ed,1),1);
            
            Edm=Ed.*Em;
            
            C=zeros(numel(Edm),1);
            [~,vSort]=sort(Edm);
            k=fix(numel(Edm)*p);
            %C(vSort(end-k+1:end))=1;
            C(vSort(end-k+1:end))=Edm(vSort(end-k+1:end));
            
            Z=C.*z;
            Z=reshape(Z,size(Km,1),size(Kd,1));
         
            matPredict=Ud*Z'*Um';
         end
         
         if Option==2
            z=reshape(Um'*matOrig'*Ud,numel(Um'*matOrig'*Ud),1);
         
            Ed=diag(ed);
            Ed=repmat(Ed',size(em,1),1);
            Ed=reshape(Ed,numel(Ed),1);
            Em=diag(em);
            Em=repmat(Em,size(ed,1),1);
            
            Edm=Ed+Em;
            
            C=zeros(numel(Edm),1);
            [~,vSort]=sort(Edm);
            k=fix(numel(Edm)*p);
            %C(vSort(end-k+1:end))=1;
            C(vSort(end-k+1:end))=Edm(vSort(end-k+1:end));
            
            Z=C.*z;
            Z=reshape(Z,size(Km,1),size(Kd,1));
         
            matPredict=Ud*Z'*Um';
         end
         
         if Option==3
            [~,vSortD]=sort(ed);
            kd=fix(size(ed,1)*p);
            Cd=zeros(size(ed,1),1);
            %Cd(vSortD(end-kd+1:end))=1;
            Cd(vSortD(end-kd+1:end))=ed(vSortD(end-kd+1:end));
            matPredictD=Ud*diag(Cd)*Ud'*matOrig;
            
            [~,vSortM]=sort(em);
            km=fix(size(em,1)*p);
            Cm=zeros(size(em,1),1);
            %Cm(vSortM(end-km+1:end))=1;
            Cm(vSortM(end-km+1:end))=em(vSortM(end-km+1:end));
            matPredictM=Um*diag(Cm)*Um'*matOrig';
            
            matPredict=(matPredictD+matPredictM')/2;
         end
            
            