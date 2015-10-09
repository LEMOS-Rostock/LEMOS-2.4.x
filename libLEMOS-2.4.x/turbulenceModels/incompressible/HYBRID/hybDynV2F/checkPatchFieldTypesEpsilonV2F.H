if (!isType<zeroGradientFvPatchScalarField>(epsilon_.boundaryField()[patchi]))
{
    FatalErrorIn("boundary-condition evaluation")
        << epsilon_.boundaryField()[patchi].type()
        << " is the wrong epsilon patchField type for boundary-condition on patch "
        << curPatch.name() << nl
        << "    should be zeroGradient"
        << exit(FatalError);
}
