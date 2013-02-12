use probabilistic_annotationImpl;

use probabilistic_annotationServer;



my @dispatch;

{
    my $obj = probabilistic_annotationImpl->new;
    push(@dispatch, 'ProbabilisticAnnotation' => $obj);
}


my $server = probabilistic_annotationServer->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler;
