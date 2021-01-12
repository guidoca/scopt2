%% Validates atmosphere model by comparison

AtmosphereConst     = atmosphereModel('const');
AtmosphereExp       = atmosphereModel('exp');
AtmosphereTabulated = atmosphereModel('tabulated');
AtmosphereConstUs76 = atmosphereModel('us76');

fig1 = figure;
fig2 = figure;

AtmosphereConst.plot(fig1); hold on
AtmosphereExp.plot(fig1); hold on
AtmosphereTabulated.plot(fig1); hold on
AtmosphereConstUs76.plot(fig1); hold on
subplot 221
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');
subplot 222
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');
subplot 223
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');
subplot 224
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');


AtmosphereConst.plotDerivatives(fig2); hold on
AtmosphereExp.plotDerivatives(fig2); hold on
AtmosphereTabulated.plotDerivatives(fig2); hold on
AtmosphereConstUs76.plotDerivatives(fig2); hold on
subplot 221
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');
subplot 222
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');
subplot 223
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');
subplot 224
legend({AtmosphereConst.model,AtmosphereExp.model,AtmosphereTabulated.model,AtmosphereConstUs76.model},'location','best');

