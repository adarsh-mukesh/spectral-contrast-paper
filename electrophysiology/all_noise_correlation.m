%% bootstrapped data prepration
load('C:\Users\amuke\Documents\MATLAB\rss-ssa\all_rstr.mat')
its{1,1}=1:13; its{1,2}=14:26;
its{2,1}=28:31; its{2,2}=32:35;
its{3,1}=37:61; its{3,2}=63:87;
a1=all_spks;
a2=zeros(size(a1,1),15,176);
for tts=1:15
    a2(:,tts,:)=nansum(a1(:,500+(tts-1)*240+1:500+tts*240,:),2);
end
nbs=1000;
all_c1=zeros(size(a2,3),nbs,3);
all_c2=zeros(size(a2,3),nbs,3);
all_m1=cell(1,4);
all_m2=cell(1,4);
for ii=1:size(a2,3)
    for jj=1:3
        for tt=1:nbs
            tr=randi([1 length(its{jj,1})],length(its{jj,1}),1);
            
            ahs1=nanmean(a2(its{jj,1}(tr),[1:7],ii),2);
            als1=nanmean(a2(its{jj,2}(tr),[1:7],ii),2);
            ahs2=nanmean(a2(its{jj,1}(tr),[2:7],ii),2);
            als2=nanmean(a2(its{jj,2}(tr),[2:7],ii),2);
            
            ahd=a2(its{jj,2}(tr),8,ii);
            ald=a2(its{jj,1}(tr),8,ii);
            
            csi1=(nanmean([ahd;ald])-nanmean([ahs1;als1]))./(nanmean([ahd;ald])+nanmean([ahs1;als1]));
            cci1=(nanmean([ald;als1])-nanmean([ahd;ahs1]))./(nanmean([ald;als1])+nanmean([ahd;ahs1]));
            dsi1_1=(nanmean(ald)-nanmean(ahs1))./(nanmean(ald)+nanmean(ahs1));
            dsi1_2=(nanmean(ahd)-nanmean(als1))./(nanmean(ahd)+nanmean(als1));
            
            csi2=(nanmean([ahd;ald])-nanmean([ahs2;als2]))./(nanmean([ahd;ald])+nanmean([ahs2;als2]));
            cci2=(nanmean([ald;als2])-nanmean([ahd;ahs2]))./(nanmean([ald;als2])+nanmean([ahd;ahs2]));
            dsi2_1=(nanmean(ald)-nanmean(ahs2))./(nanmean(ald)+nanmean(ahs2));
            dsi2_2=(nanmean(ahd)-nanmean(als2))./(nanmean(ahd)+nanmean(als2));
            
            nmb=sign([csi1 cci1]);
            nmb(nmb==0)=1;
            nmb(nmb==-1)=0;
            op=nmb*[1;2];
            switch op
                case 0
                    all_c1(ii,tt,jj)=3;
                case 1
                    all_c1(ii,tt,jj)=4;
                case 2
                    all_c1(ii,tt,jj)=2;
                case 3
                    all_c1(ii,tt,jj)=1;
            end
            all_m1{1,1}(ii,tt,jj)=csi1;
            all_m1{1,2}(ii,tt,jj)=cci1;
            all_m1{1,3}(ii,tt,jj)=dsi1_1;
            all_m1{1,4}(ii,tt,jj)=dsi1_2;
            
            nmb=sign([csi2 cci2]);
            nmb(nmb==0)=1;
            nmb(nmb==-1)=0;
            op=nmb*[1;2];
            switch op
                case 0
                    all_c2(ii,tt,jj)=3;
                case 1
                    all_c2(ii,tt,jj)=4;
                case 2
                    all_c2(ii,tt,jj)=2;
                case 3
                    all_c2(ii,tt,jj)=1;
            end
            all_m2{1,1}(ii,tt,jj)=csi2;
            all_m2{1,2}(ii,tt,jj)=cci2;
            all_m2{1,3}(ii,tt,jj)=dsi2_1;
            all_m2{1,4}(ii,tt,jj)=dsi2_2;
        end
    end
end

kts=16.*[0:11];
a1=all_spks;
a3=zeros(size(a1,1),900,176);
for tts=1:900
    a3(:,tts,:)=nansum(a1(:,(tts-1)*5+1:tts*5,:),2);
end
kvs=nchoosek([1:16],2);
all_ncr=cell(3,1);
for ii=1:length(kts)-1;
    ii
    a2=a3(:,:,kts(ii)+1:kts(ii+1));
    for jj=1:3
        aa=a2([its{jj,1} its{jj,2}],:,:);
        for kk=1:size(kvs,1)
            ad1=aa(:,101:820,kvs(kk,1));
            ad2=aa(:,101:820,kvs(kk,2));
            am1=ad1-repmat(nanmean(ad1),size(ad1,1),1);
            am2=ad2-repmat(nanmean(ad2),size(ad2,1),1);
            c1=zeros(1,size(am1,1));
            parfor ll=1:size(am1,1)
                cc=corrcoef(am1(ll,:),am2(ll,:));
                c1(1,ll)=cc(1,2);
            end
            all_ncr{jj,1}=[all_ncr{jj,1}; [kts(ii)+kvs(kk,:) nanmean(c1)]];
        end
    end
end

%% within category NC - Fig. 3B, 3F, 4F and supplementary figure C,H,M
for op=1:3
    ty=all_c1(:,:,op);
    all_bs1=[];
    all_rc1=[];
    uc=[[1 2];[3 4];[1 4];[2 3]];
    for nn=1:nbs
        nrs=cell(1,4);
        for ii=1:4
            fn=find(ty(:,nn)==ii);
            for kk=1:length(kts)-1
                ab=fn(fn<=kts(kk+1) & fn>kts(kk));
                if length(ab)<2
                    continue
                end
                bb=nchoosek(ab,2);
                [tf,index]=ismember(all_ncr{op,1}(:,[1 2]),bb,'rows');
                nrs{1,ii}=[nrs{1,ii} all_ncr{op,1}(find(tf==1),3)'];
            end
        end
        rt=cellfun(@(w) length(w),nrs);
        all_bs1=[all_bs1;cellfun(@(w) nanmean(w),nrs)];
        all_rc1=[all_rc1; nanmean(cell2mat(nrs))];
    end

    ty=all_c2(:,:,op);
    all_bs2=[];
    all_rc2=[];
    uc=[[1 2];[3 4];[1 4];[2 3]];
    for nn=1:nbs
        nrs=cell(1,4);
        for ii=1:4
            fn=find(ty(:,nn)==ii);
            for kk=1:length(kts)-1
                ab=fn(fn<=kts(kk+1) & fn>kts(kk));
                if length(ab)<2
                    continue
                end
                bb=nchoosek(ab,2);
                [tf,index]=ismember(all_ncr{op,1}(:,[1 2]),bb,'rows');
                nrs{1,ii}=[nrs{1,ii} all_ncr{op,1}(find(tf==1),3)'];
            end
        end
        rt=cellfun(@(w) length(w),nrs);
        all_bs2=[all_bs2;cellfun(@(w) nanmean(w),nrs)];
        all_rc2=[all_rc2; nanmean(cell2mat(nrs))];
    end

    ci=0.1;
    all_bs=sort(all_bs1);
    all_rc=sort(all_rc1);
    mn_bs=nanmean(all_bs);
    mn_rc=nanmean(all_rc);
    lm_bs=zeros(2,4);
    for ii=1:4
        aa=all_bs(:,ii);
        aa=aa(~isnan(aa));
        lm_bs(1,ii)=aa(fix(nbs*ci));
        lm_bs(2,ii)=aa(nbs-fix(nbs*ci));
    end
    aa=all_rc;
    aa=aa(~isnan(aa));
    lm_rc=[aa(fix(nbs*ci));aa(nbs-fix(nbs*ci))];
    figure
    shadedErrorBar([1:4],mn_rc.*ones(1,4),repmat([(lm_rc(2)-mn_rc);(mn_rc-lm_rc(1))],1,4),'k',0.5);
    hold on
    errorbar([1:4],mn_bs,mn_bs-lm_bs(1,:),lm_bs(2,:)-mn_bs,'r');
    hold on

    all_bs=sort(all_bs2);
    all_rc=sort(all_rc2);
    mn_bs=nanmean(all_bs);
    mn_rc=nanmean(all_rc);
    lm_bs=zeros(2,4);
    for ii=1:4
        aa=all_bs(:,ii);
        aa=aa(~isnan(aa));
        lm_bs(1,ii)=aa(fix(nbs*ci));
        lm_bs(2,ii)=aa(nbs-fix(nbs*ci));
    end
    aa=all_rc;
    aa=aa(~isnan(aa));
    lm_rc=[aa(fix(nbs*ci));aa(nbs-fix(nbs*ci))];
    shadedErrorBar([1:4]+0.2,mn_rc.*ones(1,4),repmat([(lm_rc(2)-mn_rc);(mn_rc-lm_rc(1))],1,4),'c',0.5);
    hold on
    errorbar([1:4]+0.2,mn_bs,mn_bs-lm_bs(1,:),lm_bs(2,:)-mn_bs,'b');

end
%% Across category NC - Fig. 3D, 3H, 4J and supplementary E,J,Q
for op=1:3
    ty=all_c1(:,:,op);
    all_bs1=[];
    all_rc1=[];
    uc=[[1 2];[3 4];[1 4];[2 3]];
    for nn=1:nbs
        nrs=cell(1,4);
        for ii=1:4
            fn1=find(ty(:,nn)==ii);
            fn2=intersect(find(ty(:,nn)~=ii),find(ty(:,nn)>0));
            for kk=1:length(kts)-1
                ab1=fn1(fn1<=kts(kk+1) & fn1>kts(kk));
                ab2=fn2(fn2<=kts(kk+1) & fn2>kts(kk));
                if isempty(ab1) || isempty(ab2)
                    continue
                end
                bb=(combvec(ab1',ab2'))';
                iu=find(bb(:,1)>bb(:,2));
                bb(iu,:)=[bb(iu,2) bb(iu,1)];
                [tf,index]=ismember(all_ncr{op,1}(:,[1 2]),bb,'rows');
                nrs{1,ii}=[nrs{1,ii} all_ncr{op,1}(find(tf==1),3)'];
            end
        end
        rt=cellfun(@(w) length(w),nrs);
        all_bs1=[all_bs1;cellfun(@(w) nanmean(w),nrs)];
        all_rc1=[all_rc1; nanmean(cell2mat(nrs))];
    end
    
    ty=all_c2(:,:,op);
    all_bs2=[];
    all_rc2=[];
    uc=[[1 2];[3 4];[1 4];[2 3]];
    for nn=1:nbs
        nrs=cell(1,4);
        for ii=1:4
            fn1=find(ty(:,nn)==ii);
            fn2=intersect(find(ty(:,nn)~=ii),find(ty(:,nn)>0));
            for kk=1:length(kts)-1
                ab1=fn1(fn1<=kts(kk+1) & fn1>kts(kk));
                ab2=fn2(fn2<=kts(kk+1) & fn2>kts(kk));
                if isempty(ab1) || isempty(ab2)
                    continue
                end
                bb=(combvec(ab1',ab2'))';
                iu=find(bb(:,1)>bb(:,2));
                bb(iu,:)=[bb(iu,2) bb(iu,1)];
                [tf,index]=ismember(all_ncr{op,1}(:,[1 2]),bb,'rows');
                nrs{1,ii}=[nrs{1,ii} all_ncr{op,1}(find(tf==1),3)'];
            end
        end
        rt=cellfun(@(w) length(w),nrs);
        all_bs2=[all_bs2;cellfun(@(w) nanmean(w),nrs)];
        all_rc2=[all_rc2; nanmean(cell2mat(nrs))];
    end
    
    ci=0.1;
    all_bs=sort(all_bs1);
    all_rc=sort(all_rc1);
    mn_bs=nanmean(all_bs);
    mn_rc=nanmean(all_rc);
    lm_bs=zeros(2,4);
    for ii=1:4
        aa=all_bs(:,ii);
        aa=aa(~isnan(aa));
        lm_bs(1,ii)=aa(fix(nbs*ci));
        lm_bs(2,ii)=aa(nbs-fix(nbs*ci));
    end
    aa=all_rc;
    aa=aa(~isnan(aa));
    lm_rc=[aa(fix(nbs*ci));aa(nbs-fix(nbs*ci))];
    figure
    shadedErrorBar([1:4],mn_rc.*ones(1,4),repmat([(lm_rc(2)-mn_rc);(mn_rc-lm_rc(1))],1,4),'k',0.5);
    hold on
    errorbar([1:4],mn_bs,mn_bs-lm_bs(1,:),lm_bs(2,:)-mn_bs,'r');
    hold on
    
    all_bs=sort(all_bs2);
    all_rc=sort(all_rc2);
    mn_bs=nanmean(all_bs);
    mn_rc=nanmean(all_rc);
    lm_bs=zeros(2,4);
    for ii=1:4
        aa=all_bs(:,ii);
        aa=aa(~isnan(aa));
        lm_bs(1,ii)=aa(fix(nbs*ci));
        lm_bs(2,ii)=aa(nbs-fix(nbs*ci));
    end
    aa=all_rc;
    aa=aa(~isnan(aa));
    lm_rc=[aa(fix(nbs*ci));aa(nbs-fix(nbs*ci))];
    shadedErrorBar([1:4]+0.2,mn_rc.*ones(1,4),repmat([(lm_rc(2)-mn_rc);(mn_rc-lm_rc(1))],1,4),'c',0.5);
    hold on
    errorbar([1:4]+0.2,mn_bs,mn_bs-lm_bs(1,:),lm_bs(2,:)-mn_bs,'b');
end
%% CCN, CCX NC - Fig. 3C, 3G, 4I and supplementary D,I,H
for op=1:3
    ty=all_c1(:,:,op);
    all_bs1=[];
    all_rc1=[];
    uc=[[1 2];[3 4];[1 4];[2 3]];
    for nn=1:nbs
        nrs=cell(1,4);
        for ii=1:4
            %fn=find(ty(:,nn)==ii);
            fn=union(find(ty(:,nn)==uc(ii,1)),find(ty(:,nn)==uc(ii,2)));
            for kk=1:length(kts)-1
                ab=fn(fn<=kts(kk+1) & fn>kts(kk));
                if length(ab)<2
                    continue
                end
                bb=nchoosek(ab,2);
                [tf,index]=ismember(all_ncr{op,1}(:,[1 2]),bb,'rows');
                nrs{1,ii}=[nrs{1,ii} all_ncr{op,1}(find(tf==1),3)'];
            end
        end
        rt=cellfun(@(w) length(w),nrs);
        all_bs1=[all_bs1;cellfun(@(w) nanmean(w),nrs)];
        all_rc1=[all_rc1; [nanmean(cell2mat(nrs(1:2))) nanmean(cell2mat(nrs(3:4)))]];
    end
    
    ty=all_c2(:,:,op);
    all_bs2=[];
    all_rc2=[];
    uc=[[1 2];[3 4];[1 4];[2 3]];
    for nn=1:nbs
        nrs=cell(1,4);
        for ii=1:4
            %fn=find(ty(:,nn)==ii);
            fn=union(find(ty(:,nn)==uc(ii,1)),find(ty(:,nn)==uc(ii,2)));
            for kk=1:length(kts)-1
                ab=fn(fn<=kts(kk+1) & fn>kts(kk));
                if length(ab)<2
                    continue
                end
                bb=nchoosek(ab,2);
                [tf,index]=ismember(all_ncr{op,1}(:,[1 2]),bb,'rows');
                nrs{1,ii}=[nrs{1,ii} all_ncr{op,1}(find(tf==1),3)'];
            end
        end
        rt=cellfun(@(w) length(w),nrs);
        all_bs2=[all_bs2;cellfun(@(w) nanmean(w),nrs)];
        all_rc2=[all_rc2; [nanmean(cell2mat(nrs(1:2))) nanmean(cell2mat(nrs(3:4)))]];
    end
    
    ci=0.1;
    all_bs=sort(all_bs1);
    all_rc=sort(all_rc1);
    mn_bs=nanmean(all_bs);
    mn_rc=nanmean(all_rc);
    lm_bs=zeros(2,4);
    for ii=1:4
        aa=all_bs(:,ii);
        aa=aa(~isnan(aa));
        lm_bs(1,ii)=aa(fix(nbs*ci));
        lm_bs(2,ii)=aa(nbs-fix(nbs*ci));
    end
    lm_rc=[];
    for ii=1:2
        aa=all_rc(:,ii);
        aa=aa(~isnan(aa));
        lm_rc=[lm_rc [aa(fix(nbs*ci));aa(nbs-fix(nbs*ci))]];
    end
    figure
    shadedErrorBar([1:2],mn_rc(1)*ones(1,2),repmat([(lm_rc(2,1)-mn_rc(1));(mn_rc(1)-lm_rc(1,1))],1,2),'k',0.5);
    hold on
    shadedErrorBar([3:4],mn_rc(2)*ones(1,2),repmat([(lm_rc(2,2)-mn_rc(2));(mn_rc(2)-lm_rc(1,2))],1,2),'k',0.5);
    hold on
    errorbar([1:4],mn_bs,mn_bs-lm_bs(1,:),lm_bs(2,:)-mn_bs,'r');
    hold on
    
    ci=0.1;
    all_bs=sort(all_bs2);
    all_rc=sort(all_rc2);
    mn_bs=nanmean(all_bs);
    mn_rc=nanmean(all_rc);
    lm_bs=zeros(2,4);
    for ii=1:4
        aa=all_bs(:,ii);
        aa=aa(~isnan(aa));
        lm_bs(1,ii)=aa(fix(nbs*ci));
        lm_bs(2,ii)=aa(nbs-fix(nbs*ci));
    end
    lm_rc=[];
    for ii=1:2
        aa=all_rc(:,ii);
        aa=aa(~isnan(aa));
        lm_rc=[lm_rc [aa(fix(nbs*ci));aa(nbs-fix(nbs*ci))]];
    end
    shadedErrorBar([1:2]+0.2,mn_rc(1)*ones(1,2),repmat([(lm_rc(2,1)-mn_rc(1));(mn_rc(1)-lm_rc(1,1))],1,2),'c',0.5);
    hold on
    shadedErrorBar([3:4]+0.2,mn_rc(2)*ones(1,2),repmat([(lm_rc(2,2)-mn_rc(2));(mn_rc(2)-lm_rc(1,2))],1,2),'c',0.5);
    hold on
    errorbar([1:4]+0.2,mn_bs,mn_bs-lm_bs(1,:),lm_bs(2,:)-mn_bs,'b');
end
%% Fig 3i, 3j and supplementary R,S
nbs=1000;
nrs1=zeros(1000,4);
nrs2=zeros(1000,4);
bss=zeros(1000,4);
for ii=1:1000
    a1=squeeze(all_c1(:,ii,1:2)); % all_c1 for Fig.3i,3j , all_c2 for supplementary R, S
    cl1=find(a1(:,2)==4);% == 3 for 3i and supplementary R, == 4 for 3j and supplementary S
    cts1=[cl1 a1(cl1,1)];
    for jj=1:4
        kbs=cts1(find(cts1(:,2)==jj),1);
        [tf1,index]=ismember(all_ncr{1,1}(:,1),kbs);
        [tf2,index]=ismember(all_ncr{1,1}(:,2),kbs);
        nrs2(ii,jj)=nanmean(all_ncr{1,1}(union(find(tf1==1),find(tf2==1)),3));
        [tf3,index]=ismember(all_ncr{2,1}(:,1),kbs);
        [tf4,index]=ismember(all_ncr{2,1}(:,2),kbs);
        nrs1(ii,jj)=nanmean(all_ncr{2,1}(union(find(tf3==1),find(tf4==1)),3));
        bss(ii,jj)=nanmean([all_ncr{1,1}(union(find(tf1==1),find(tf2==1)),3);all_ncr{2,1}(union(find(tf3==1),find(tf4==1)),3)]);
    end
end
ci=0.1;
nrs1=sort(nrs1);
mn1=nanmean(nrs1);
ll1=nrs1(fix(ci*nbs),:);
ul1=nrs1(nbs-fix(ci*nbs),:);
nrs2=sort(nrs2);
mn2=nanmean(nrs2);
ll2=nrs2(fix(ci*nbs),:);
ul2=nrs2(nbs-fix(ci*nbs),:);
bss=sort(bss);
mbs=nanmean(bss);
llb=bss(fix(ci*nbs),:);
ulb=bss(nbs-fix(ci*nbs),:);
figure
for ii=1:4
    shadedErrorBar([ii ii+0.2],ones(1,2)*mbs(ii),[ones(1,2)*(ulb(ii)-mbs(ii));ones(1,2)*(mbs(ii)-llb(ii))],'k');
    hold on
end
errorbar([1:4],mn1,mn1-ll1,ul1-mn1,'b');
hold on
errorbar([1:4]+0.2,mn2,mn2-ll2,ul2-mn2,'r');
hold on










