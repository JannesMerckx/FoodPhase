%% setup

n = 12;
[X,Y] = meshgrid(linspace(0,1,n),linspace(0,1,n));
x = [X(:),Y(:)];

%refractive index
eta = ones(n_2,1) + .3*rand(n_2,1);

%parameters
PLOT = true;
RAY_TEST = false;



%% triangulation and mesh admninistration
n_0 = size(x,1);
DT = delaunay(x);
n_2 = size(DT,1);

if PLOT
    close all
    triplot(DT,x(:,1),x(:,2));
    axis square
end


cEL = zeros(n_2,3);
vCL = zeros(n_0,6);
EL = zeros(3*n_2,2);
n_1 = 0;
for cell = 1:n_2
    
    for l=1:3
        
        edge = [DT(cell,l),DT(cell,mod(l,3)+1)];
        
        inlist = false;
        f = 0;
        for jj=1:n_1
            if norm(EL(jj,:) - edge)==0
                inlist = true;
                break
            elseif norm(EL(jj,:) - fliplr(edge))==0
                inlist = true;
                break
            end
        end
        
        if inlist == true
            cEL(cell,l) = jj;
            continue
        end
        
        
        
        n_1 = n_1+1;
        EL(n_1,:) = edge;
        cEL(cell,l) = n_1;
    end
    
    
    verts = DT(cell,:);
    for i=1:3
        for j=1:6
            if vCL(verts(i),j) == 0
                vCL(verts(i),j) = cell;
                break
            end
        end
    end
    
    
end
EL = EL(1:n_1,:);





B21 = sparse(n_2,n_1);
B10 = sparse(n_1,n_0);

for cell = 1:n_2
    e = cEL(cell,:);
    v = DT(cell,:);
    
    vec1 = x(EL(e(1),2),:)-x(EL(e(1),1),:);
    vec2 = x(EL(e(2),2),:)-x(EL(e(2),1),:);
    vec3 = x(EL(e(3),2),:)-x(EL(e(3),1),:);
    c = (x(DT(cell,1),:)+x(DT(cell,2),:)+x(DT(cell,3),:))/3;
    
    for i=1:3
        B10(e(i),EL(e(i),1))=-1;
        B10(e(i),EL(e(i),2))=1;
    end
    
    
    
    %f = [c(2)-y;x-c(1)];
    
    f1 = [c(2)-x(EL(e(1),1),2);x(EL(e(1),1),1)-c(1)];
    f2 = [c(2)-x(EL(e(2),1),2);x(EL(e(2),1),1)-c(1)];
    f3 = [c(2)-x(EL(e(3),1),2);x(EL(e(3),1),1)-c(1)];
    
    
    
    
    s = sign([vec1*f1,vec2*f2,vec3*f3]);
    if any(s==0)
        error('test vector field error');
    end
    B21(cell,e) = s;
    
    
end



%% test ray
if RAY_TEST
    th = 0;
    delta = .0001;
    
    
    for ray = 0:1
        
        r = [-sin(th),cos(th),-.5*cos(th)+ray*delta];
        f = 1;
        
        
        
        int = [0,-r(3)/r(2)];
        
        
        
        hold on
        plot(int(1),int(2),'kx','MarkerSize',8)
        
        
        
        
        [~,j] = min( (x(:,1)-int(1)).^2 + (x(:,2)-int(2)).^2  );
        
        
        
        edges = find(B10(:,j))';
        
        for e = edges
            v1 = int;
            v2 = x(EL(e(1),1),:);
            v3 = x(EL(e(1),2),:);
            s = SA(v1,v2,v3);
            
            
            on_edge = false;
            if abs(s)<1e-13
                lambda = (v2-v1)./(v2-v3);
                if isnan(lambda(1))
                    lambda = lambda(2);
                else
                    lambda = lambda(1);
                end
                
                if 0<= lambda && lambda <= 1
                    on_edge = true;
                    break
                end
                
            end
        end
        
        
        hold on
        %plot([x(EL(e,1),1),x(EL(e,2),1)],[x(EL(e,1),2),x(EL(e,2),2)],'g-');
        
        
        
        ce = e;
        cc = find(B21(:,e));
        
        
        
        p = eta(cc)*[r(2),-r(1)];
        
        
        
        while max(abs(int))<1
            
            %search edges
            edges = setdiff(cEL(cc,:),ce);
            
            for edge = edges
                
                x1 = x(EL(edge,1),:);
                x2 = x(EL(edge,2),:);
                t = x2-x1;
                t = t/norm(t);
                n = [-t(2),t(1)];
                le = [n(1),n(2),n*x2'];
                
                if abs(r(1)*le(2)-r(2)*le(1))<1e-13
                    %ray and edge are parrallel
                    continue
                end
                
                %compute new intersection point parameters
                %ray in line form
                %edge in parameter form
                lambda = r(1:2)*(int-x1)'/(r(1:2)*(x2-x1)');
                intn = (1-lambda)*x1 + lambda*x2;
                
                
                
                if (0 <= lambda && lambda <= 1)
                    %ray hits edge line segment
                    break
                end
                
                
            end
            
            %compute new cell
            cells = find(B21(:,edge))';
            ccn = setdiff(cells,cc);
            
            
            %plot([x2(1),x1(1)],[x2(2),x1(2)],'g-','MarkerSize',8);
            %plot(intn(1),intn(2),'gx','MarkerSize',8);
            plot([intn(1),int(1)],[intn(2),int(2)],'g-');
            int = intn;
            
            if isempty(ccn)
                break
            end
            
            
            %pgon = polyshape(x(DT(ccn,:),1)',x(DT(ccn,:),2)');
            %plot(pgon);
            
            
            
            
            
            pt = (p*t');
            d = eta(ccn)^2-pt^2;
            if d>=0
                pn = pt*t + sign(p*n')*sqrt(eta(ccn)^2-pt^2)*n;
            else
                pn = p - 2*(p*n')*n;
                ccn = cc;
            end
            
            
            f = f*eta(cc)/eta(ccn)*(pn*n')/(p*n');
            
            
            
            
            
            
            p = pn;
            cc = ccn;
            ce = edge;
            r = [-p(2),p(1),-p(2)*int(1)+p(1)*int(2)];
            
            
            
        end
        
        f_end = f;
        f = -f*(p*n')/eta(cc);
        
        plot(int(1),int(2),'kx','MarkerSize',8)
        axis equal
        
    end
end


%% Eikonal solver


hold on
S = inf(n_0,1);
current = logical(x(:,1)<1e-13);
%current = sparse(false(n_0,1));
%current(randi(n_0)) = 1;
S(current) = 0;
unvisited = ~current;
visited = sparse(false(n_0,1));
S(current) = 0;



if sum(current)==1 %special case
    source = find(current);
    neigh = any(B10(logical(B10(:,source)),:),1) & (~current');
    
    
    for j = find(neigh)
        cells = find(B21(:,B10(:,source) & B10(:,j)));
        ri = eta(cells(1));
        
        S(j) = norm(x(j,:)-x(source,:))*ri;
        
        current(j) = true;
    end
    
    current(source) = false;
    visited(source) = true;
end





while any(unvisited)
    
    %find vertex closest to front
    current_n = current;
    edge = any(B10(:,current),2);
    neigh = any(B10(edge,:),1)' & ~current & ~ visited;
    front = any(B10(:,current),2) & ~any(B10(:,neigh),2) & ~any(B10(:,visited),2);
    cfront = any(B21(:,front),2);
    
    [~,j] = min(S(current));
    v = find(current);
    u = v(j);
    
    
    
    ce = B10(:,u) & front;
    
    for c=find(any(B21(:,ce),2))'
        
        vert = any(B10(any(B21(c,:),1)',:),1)';
        if any(vert & visited)
            continue
        end
        w = setdiff(find(vert & current),u);
        v = setdiff(find(vert),[u;w]);
        
        if PLOT
            hold off
            triplot(DT,x(:,1),x(:,2));
            hold on
            axis([0 1 0 1])
            axis square
            
            
            %     for j=find(neigh)
            %         plot(x(j,1),x(j,2),'rx')
            %     end
            
            for j=find(current)'
                plot(x(j,1),x(j,2),'gx')
            end
            
            
            
            
            %     for e = find(edge)'
            %         plot(x(EL(e,:),1),x(EL(e,:),2),'g-');
            %     end
            
            for e = find(front)'
                plot(x(EL(e,:),1),x(EL(e,:),2),'r-');
            end
            
            pol = polyshape(x(DT(c,:),1),x(DT(c,:),2));
            plot(pol)
            axis([0 1 0 1])
            axis square
            
            
            plot(x(v,1),x(v,2),'r^')
            plot(x(w,1),x(w,2),'rs')
            plot(x(u,1),x(u,2),'ro')
        end
        
        if isempty(v) || any(isinf(S([u;w])))
            continue
        end
        
        
        
        S_lim = min(S(u)+eta(c)*norm(x(u,:)-x(v,:)),S(w)+eta(c)*norm(x(w,:)-x(v,:)));
        
        xl = @(lam) (1-lam)*x(u,:) + lam*x(w,:);
        f = @(lam) (1-lam)*S(u) + lam*S(w) + eta(c)*norm(xl(lam)-x(v,:));
        
        df = @(lam) S(w)-S(u) + eta(c)./norm(xl(lam)-x(v,:))...
            .*((x(w,:)-x(u,:))*(xl(lam)-x(v,:))');
        
        
        test = df(0)*df(1);
        if test>0
            if f(0) < f(1)
                S(v) = f(0);
                xs = xl(0);
            else
                S(v) = f(1);
                xs = xl(1);
            end
            
            if PLOT
                plot(xs(1),xs(2),'gx');
                plot([x(v,1),xs(1)],[x(v,2),xs(2)],'g-');
                axis([0 1 0 1])
                axis square
            end
            
            
        else
            lambda = fzero(df,[0 1]);
            St = f(lambda);
            if St<S(v)
                S(v) = St;
            end
            
            
            
            xs = xl(lambda);
            if PLOT
                plot(xs(1),xs(2),'gx');
                plot([x(v,1),xs(1)],[x(v,2),xs(2)],'g-');
                axis([0 1 0 1])
                axis square
            end
            
            
            
            
            
            if S(v)>=S_lim
                
                'stahp';
            end
        end
        
        current_n(v) = true;
        drawnow()
        
        
        
        
    end
    
    
    
    if any(current) == false
        visited = true(n_0,1);
        unvisited = sparse(false(n_0,1));
    else
        for u = find(current)'
            neigh = setdiff(find(any(B10(any(B10(:,u),2),:),1)),u);
            
            if all(current(neigh) | visited(neigh))
                current_n(u) = false;
                visited(u) = true;
            end
        end
    end
    current = current_n;
    
    
    
    
    
    
end





if PLOT
hold off
trisurf(DT,x(:,1),x(:,2),S)
end

%% Backward rays

%reconstruct a ray at the right edge and propagate backwards.




%% functions

function s = SAO(v1,v2)
%computes the signed area of the triangle v1,v2,O
s = .5*(v1(1)*v2(2)-v1(2)*v2(1));
end


function s = SA(v1,v2,v3)
%computes the signed area of the triangle v1,v2,v3
s = SAO(v1,v2)+SAO(v2,v3)+SAO(v3,v1);
end




