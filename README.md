# B2a-いじった場所
・DetectorConstructionでWorldの配置
targetは磁場空間
trackerは検出器としてポリスチレン


・PrimaryGeneretorで線源の定義
Sr線源の定義
if (fParticleGun->GetParticleDefinition()==G4Geantino::Geantino()){
    G4int Z = 38;
    G4int A = 90;
    G4double ionCharge = 0*eplus;
    G4double excitEnergy = 0*keV;

    

    G4ParticleDefinition* ionDefinition 
      = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ionDefinition);
    fParticleGun->SetParticleCharge(ionCharge);
 }
 線源を定義するにはBiasedRDPhiscsとPhysicslistが必要だった

 ・TrackerHit.hhでSetとGetの定義
 ・TrackerHit.ccでほしい情報(fMomとか)をいろいろやる
コンストラクタ (B2TrackerHit::B2TrackerHit()):fTrackID、fChamberNb、fEdep、fPos、および fMom などのメンバ変数を初期化します。G4VHit を継承しています。
デストラクタ (B2TrackerHit::~B2TrackerHit()):空のデストラクタ。
コピーコンストラクタ (B2TrackerHit::B2TrackerHit(const B2TrackerHit& right)):既存の B2TrackerHit オブジェクトに基づいて新しいオブジェクトを初期化するコピーコンストラクタ。
代入演算子 (const B2TrackerHit& B2TrackerHit::operator=(const B2TrackerHit& right)):代入演算子のオーバーロード。一つの B2TrackerHit オブジェクトの値を別のオブジェクトに代入します。
等号演算子 (G4bool B2TrackerHit::operator==(const B2TrackerHit& right) const):等号演算子のオーバーロード。2つの B2TrackerHit オブジェクトを比較して等しいかどうかを判断します。
Draw() メソッド:Geant4の視覚化を使用して、ヒットの位置を表す円を描画します。
Print() メソッド:コンソールに fTrackID、fChamberNb、fEdep、fPos、および fMom に関する情報を出力します。

・TrackerSDで情報の取得：edep>=0の時のhitごとに情報の取得
・RunActionで情報を入れるコラムの作成、file.rootに入れていく
・EventActionでコラムに情報を入れていく
・vis.macでγ線とか消せる
