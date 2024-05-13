@testset "types" begin
    m = MTModel([100., 100.], [100.])
    @test typeof(m) <: AbstractGeophyModel
    @test typeof(zero(m)) <: AbstractGeophyModel
    @test typeof(copy(m)) <: AbstractGeophyModel


    ω = [1., 10.];
    resp = forward(m, ω);

    @test typeof(resp) <: AbstractGeophyResponse
    @test typeof(zero(resp)) <: AbstractGeophyResponse
    @test typeof(copy(resp)) <: AbstractGeophyResponse

end