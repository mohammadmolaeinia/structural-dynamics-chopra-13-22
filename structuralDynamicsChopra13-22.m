clc; clear; close all;

%%
% 1-1) Calclute the Mass matrix M (kg) and Calclute the Stiffness matrix K (N/m)

disp('Part 1:');

m1 = input('M1(kg)(Exmp=114750):');
m2 = input('M2(kg)(Exmp=114750):');
m3 = input('M3(kg)(Exmp=114.750):');

M = [m1 0 0;0 m2 0;0 0 m3];

disp('1) Mass matrix M (kg)'); disp(M);

k1 = input('K1(N/m)(Exmp=9000000):');
k2 = input('K2(N/m)(Exmp=9000000):');
k3 = input('K3(N/m)(Exmp=900.0000):');

K = [ 1/(3*k1) , 5/(6*k1) , 4/(3*k1) ; 5/(6*k1) , 8/(3*k1) , 14/(3*k1) ;  4/(3*k1) ,  14/(3*k1) ,   26/(3*k1)+ 1/(3*k3) ]^(-1);

disp('2) Stiffness matrix K (N/m)'); disp(K);

runAgain = true;

while runAgain

    disp('Part 3: p = [1000 ; 0 ; 500]*4448.22');
    disp('Part 4: p = [1000 ; 0 ; 500]*4448.22*sin(w1*t)');
    disp('Part 5: ElCentro Modal Analysis');
    disp('Part 6: ElCentro Spectral Analysis');
    partNum = input('put Part Number 3, 4, 5, 6: ');

    %%
    % 1-3) findig phi, Frequencies, Period Time and Omega Matrix and Ploting Mode shapes

    [phi, omega2Mat] = eig(K, M);

    omega = sqrt([omega2Mat(1,1) ; omega2Mat(2,2) ; omega2Mat(3,3)]); % rad/s

    f     = omega/(2*pi);          % Hz
    T_n   = 1./f;                  % Seconds

    %%
    % 1-4) Calclute the third xi and C matrix

    %xi1 = input('xi1(Exmp=0.05):');
    %xi2 = input('xi2(Exmp=0.05):');

    xi1 = 0.05;
    xi3 = 0.05;

    xi = [xi1; xi3];

    ab = 0.5*[1/omega(1) omega(1) ; 1/omega(3) omega(3)]\xi;
    alpha = ab(1);
    beta  = ab(2);

    C = alpha*M + beta*K;

    zeta = zeros(3,1);
    for r = 1:3
        zeta(r) = 0.5*( alpha*(1/omega(r)) + beta*omega(r) );
    end

    %%
    % 1-5) Calclute Effective Modal Mass and Height

    hVector = [3; 6; 9];
    mStar = zeros(3,1);
    hStar = zeros(3,1);
    iota = [1; 1; 1];

    for r = 1:3

        Ln = phi(:,r)' * M * iota;
        Lteta = (phi(:,r)' * M * hVector);
        Mn = phi(:,r)' * M * phi(:,r);
        mStar(r) = (Ln^2) / Mn;
        hStar(r) = Lteta / Ln;
    end

    %%
    % 1-6) Displaying 1, 2 and 7 Answers and Figures

    disp('Part 2:');
    disp('1) Damping matrix C (NÂ·s/m)'); disp(C);
    disp('2) Natural periods T (s)'); disp(T_n);
    disp('3) Natural frequencies f (Hz)'); disp(f);
    disp('4) Natural frequencies omega (rad/s)'); disp(omega);
    disp('5) Mode shapes phi'); disp(phi);
    disp('6) Modal damping ratios zeta'); disp(zeta);
    disp('Part 7:');
    disp('1) Effective Modal Masses M*n (kg)'); disp(mStar);
    disp('2) Effective Modal Heights h*n (m)'); disp(hStar);

    phiNormForFigure = [ 0 0 0 ; 1 0 0 ; 0 1 0 ; 0 0 1] * phi;
    Stories = [0 1 2 3];

    figure('Name','2) The Mode shapes','NumberTitle','off')
    ax1 = subplot(1,3,1);
    plot(ax1,phiNormForFigure(:,1),Stories)
    title(ax1,'phi 1')
    ylabel(ax1,'Stories')
    xlabel(ax1,'Displacement')

    ax2 = subplot(1,3,2);
    plot(phiNormForFigure(:,2),Stories)
    title(ax2,'phi 2')
    ylabel(ax2,'Stories')
    xlabel(ax2,'Displacement')

    ax3 = subplot(1,3,3);
    plot(phiNormForFigure(:,3),Stories)
    title(ax3,'phi 3')
    ylabel(ax3,'Stories')
    xlabel(ax3,'Displacement')

    %%
    % NEWMARK’S METHOD and LINEAR SYSTEMS for multi degree of freedom systems
    % 1) Initial calculations

    ElCentro = load('ElCentro.txt');
    uddg = (0.65 / max(ElCentro) * 9.81 * ElCentro');
    dt = 0.02;
    totalTime = dt*(length(uddg)-1);
    t  = 0:dt:totalTime;

    if partNum == 3

        p = [1000 ; 0 ; 500]*ones(1,length(t))*4448.22;

    elseif partNum == 4

        p = [1000 ; 0 ; 500]*4448.22*sin(omega(1)*t);

    elseif partNum == 5

        p = -M * iota * uddg;

    elseif partNum == 6

        p = 6;

    else

        p = 6;
        disp('False input!');

    end


    if p ~= 6

        gamma = 0.5;
        beta = 0.25;

        q   = zeros(3,length(uddg));
        qd  = zeros(3,length(uddg));
        qdd = zeros(3,length(uddg));
        q(:,1) = 0;
        qd(:,1) = 0;

        P0 = p(:,1);
        qdd(:,1) = M \ (P0 - (C * qd(:,1)) - (K * q(:,1)));

        a1 = (1/(beta*dt^2))*M+((gamma)/(beta*dt))*C;
        a2 = (M/(beta*dt))+((gamma/beta-1)*C);
        a3 = (1/(beta*2)-1)*M+(dt*((gamma*0.5/beta)-1)*C);

        KEff = K + a1;

        % 2) Calculations for each time step

        for i = 1:length(uddg)-1

            PEff = (p(:,i+1)) + a1*q(:,i) + a2*qd(:,i) + a3*qdd(:,i);

            q(:,i+1) = KEff \ PEff;

            qd(:,i+1) = (gamma/(beta*dt)*(q(:,i+1)-q(:,i)))+((1-gamma/beta)*qd(:,i))+(dt*(1-gamma*0.5/beta)*qdd(:,i));

            qdd(:,i+1) = (1/(beta*dt^2)*(q(:,i+1)-q(:,i)))-(1/(beta*dt)*qd(:,i))-((0.5/beta-1)*qdd(:,i));

        end

        %%
        % A) Calculations Displacement History and top floore's drift

        u = phi * q;
        drift01 = u(1,:);
        drift12 = u(2,:) - u(1,:);
        drift23 = u(3,:) - u(2,:);

        % B) Calculations V Base and M Base

        vBase = [1 1 1] * K * u;
        mBase = [3 6 9] * K * u;

        % C) Calculations Max Dis and Max Drift

        maxU = [0 0 0;1 0 0;0 1 0;0 0 1]*[ max(abs(u(1,:))) ; max(abs(u(2,:))) ; max(abs(u(3,:))) ];
        maxDrift = [ max(abs(drift01)) ; max(abs(drift12)) ;max(abs(drift23)) ];

        % D) Calculations Max V Base and Max M Base

        maxVBase = max(abs(vBase));
        maxMBase = max(abs(mBase));

        % A-f)Displaying 3, 4, 5 and 6 A Part Figures:

        figure('Name','A) The Story Displacement and Top Floor Drift History','NumberTitle','off')
        ax1 = subplot(4,1,1);
        plot(ax1,t,u(1,:))
        title(ax1,'S1')
        xlabel(ax1,'time')
        ylabel(ax1,'Displacement')

        ax2 = subplot(4,1,2);
        plot(ax2,t,u(2,:))
        title(ax2,'S2')
        xlabel(ax2,'time')
        ylabel(ax2,'Displacement')

        ax3 = subplot(4,1,3);
        plot(ax3,t,u(3,:))
        title(ax3,'S3')
        xlabel(ax3,'time')
        ylabel(ax3,'Displacement')

        ax4 = subplot(4,1,4);
        plot(ax4,t,drift23)
        title(ax4,'Top flor Drift')
        xlabel(ax4,'time')
        ylabel(ax4,'Drift')

        % B-f)Displaying 3, 4, 5 and 6 B Part Figures:

        figure('Name','B) The V Base and M Base History','NumberTitle','off')
        ax1 = subplot(2,1,1);
        plot(ax1,t,vBase/1000)
        title(ax1,'V Base')
        xlabel(ax1,'time')
        ylabel(ax1,'V Base kN')

        ax2 = subplot(2,1,2);
        plot(ax2,t,mBase/1000)
        title(ax2,'M Base')
        xlabel(ax2,'time')
        ylabel(ax2,'M base kN.m')

        % C-f)Displaying 3, 4 and 5 C Part Figures:

        figure('Name','C) Max Dis and Max Drift','NumberTitle','off')
        ax1 = subplot(1,2,1);
        plot(ax1,maxU,Stories)
        title(ax1,'Max Displacement of Stories')
        ylabel(ax1,'Stories')
        xlabel(ax1,'Displacement')

        ax2 = subplot(1,2,2);
        plot(ax2,maxDrift,[1 2 3])
        title(ax2,'Max Drift of floor')
        xlabel(ax2,'Drift')
        ylabel(ax2,'floor')

        % D-f)Displaying 3, 4 and 5 D Part Answer:

        disp('Part 3, 4 ,5:');
        disp('1) Max V Base (kN)'); disp(maxVBase/1000);
        disp('2) Max M Base (kN.m)'); disp(maxMBase/1000);

    else

        TnRS = 0:0.01:5;
        RS = zeros(3,length(TnRS));

        for z = 1:3

            for Tn = 1:length(TnRS)

                xiRS = zeta(z);
                wnRS = 2*pi/TnRS(Tn);
                mRS = 1;
                kRS = wnRS^2;
                cRS = 2*wnRS*xiRS^2;

                pRS = -mRS * ElCentro';

                gammaRS = 0.5;
                betaRS = 0.25;

                uRS   = zeros(1,length(ElCentro'));
                udRS  = zeros(1,length(ElCentro'));
                uddRS = zeros(1,length(ElCentro'));

                udd0RS = (pRS(1) - cRS * udRS(1) - kRS * uRS(1))/mRS;

                a1RS = (1/(betaRS*dt^2))*mRS+((gammaRS)/(betaRS*dt))*cRS;
                a2RS = (mRS/(betaRS*dt))+((gammaRS/betaRS-1)*cRS);
                a3RS = (1/(betaRS*2)-1)*mRS+(dt*((gammaRS*0.5/betaRS)-1)*cRS);

                kEffRS = kRS + a1RS;

                for i = 1:length(ElCentro')-1

                    PEffRS = (pRS(i+1)) + a1RS*uRS(i) + a2RS*udRS(i) + a3RS*uddRS(i);

                    uRS(i+1) =PEffRS / kEffRS;

                    udRS(i+1) = (gammaRS/(betaRS*dt)*(uRS(i+1)-uRS(i)))+((1-gammaRS/betaRS)*udRS(i))+(dt*(1-gammaRS*0.5/betaRS)*uddRS(i));

                    uddRS(i+1) = (1/(betaRS*dt^2)*(uRS(i+1)-uRS(i)))-(1/(betaRS*dt)*udRS(i))-((0.5/betaRS-1)*uddRS(i));

                end

                RS(z,Tn) = max(abs(uddRS));

            end

        end

        % C) Calculations Max Dis and Max Drift

        AModal = [RS(1,round(T_n(1)/0.01 , 0)+1) ; RS(2,round(T_n(2)/0.01 , 0)+1) ; RS(3,round(T_n(3)/0.01 , 0)+1)] * 0.65 / max(ElCentro) * 9.81;
        maxq = AModal ./ ( omega .^2 );
        DModal = phi .* maxq';
        DFlor = sqrt((DModal.^2)*[1 ; 1 ; 1]);
        maxU = [0 0 0;1 0 0;0 1 0;0 0 1] * DFlor;

        driftModal = [DModal(1,:) ; DModal(2,:) - DModal(1,:) ; DModal(3,:) - DModal(2,:)].^2;
        driftFlor = sqrt((driftModal.^2)*[1 ; 1 ; 1]);
        maxDrift = driftFlor;

        % D) Calculations Max V Base and Max M Base

        fsoModal = M * (phi .* AModal');
        VBaseModal = (fsoModal')*[1 ; 1 ; 1];
        maxVBase = sqrt((VBaseModal'.^2)*[1 ; 1 ; 1]);

        MBaseModal = (fsoModal')*[3 ; 6 ; 9];
        maxMBase = sqrt((MBaseModal'.^2)*[1 ; 1 ; 1]);

        % C-f)Displaying 6 Response Spectral Part Figures:

        figure('Name','6) The Response Spectral','NumberTitle','off')
        ax1 = subplot(1,1,1);
        plot(ax1,TnRS,RS)
        title(ax1,'Response Spectral for Xi 1 and 2')
        xlabel(ax1,'Tn')
        ylabel(ax1,'A/udd')

        % C-f)Displaying 6 C Part Figures:

        figure('Name','C) Max Dis and Max Drift','NumberTitle','off')
        ax1 = subplot(1,2,1);
        plot(ax1,maxU,Stories)
        title(ax1,'Max Displacement of Stories')
        ylabel(ax1,'Stories')
        xlabel(ax1,'Displacement')

        ax2 = subplot(1,2,2);
        plot(ax2,maxDrift,[1 2 3])
        title(ax2,'Max Drift of floor')
        xlabel(ax2,'Drift')
        ylabel(ax2,'floor')

        % D-f)Displaying 6 D Part Answer:
        disp('Part 6:');

        disp('1) Max V Base (kN)'); disp(maxVBase/1000);
        disp('2) Max M Base (kN.m)'); disp(maxMBase/1000);

    end
    
    response = input('to Continue Type anything but 0 and to End 0: ');
    
    if response == 0
        runAgain = false;
    end
    
end